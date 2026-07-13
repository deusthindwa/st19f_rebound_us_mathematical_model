// by Deus Thindwa, Dan Weinberger, Ginny Pitzer
// age-structured mathematical model for pneumococcal transmission
// 05/23/2026

// =============================================================================

// Stan model for fitting PCV7 carriage-blocking effectiveness (delta_V, delta_F)
// and waning rates (omega_V, omega_F) by HMC (NUTS).
//
// demographic dynamics are  driven INSIDE Stan by `birth_rate_year[n_year]`, `mu_age[n_age]` and a single `init_pop_1999` vector. 
// the year-boundary update applies births, deaths and aging via a 12-substep Euler integration of the linear demographic ODE.
// FOI normalisation uses the per-age compartment sum at year start, recomputed each year
// the contact matrix is symmetrised inside Stan from the raw daily contacts and the current model populations. 
// pre-PCV in-Stan burn-in keeps "rescale to init_pop_1999" behaviour to bring carriage state to its pre-vaccination equilibrium at 1999 pop levels 
//
// priors:
//   delta_V, delta_F ~ Beta(13, 7)   (mean 0.70, 95 % CrI ~ 0.44-0.83)
//   omega_V, omega_F ~ Gamma(4, 11)   (mean 0.25, 95 % CrI ~ 0.030-0.697, duration mean ~4 yr)
//   log_q_age[a]  ~ Uniform(log_q_lower[a], log_q_upper[a])
//   log_rr_V      ~ Uniform(log_rr_V_lb, log_rr_V_ub)
//   log_rr_F      ~ Uniform(log_rr_F_lb, log_rr_F_ub)
//   log_rr_N      ~ Uniform(log_rr_N_lb, log_rr_N_ub)

// bounds = marginal 99th quantiles of carriage-fit posterior, samples are independent across the seven, log-uniform on the natural scale
//
// vaccination at birth: per-boundary coverage `vacc_cov_year[t]` is the
// fraction of newborns entering S_7<1y rather than S<1y. Current schedule is (0.43, 0.99, 0.99, 0.99, 0.95, ..., 0.95).
//
// state layout per age group (n_age = 4):
//   unvaccinated:  S, V, F, N, VF, NV, NF (families 0..6)
//   vaccinated:    S_7, V_7, F_7, N_7, VF_7, NV_7, NF_7 (families 7..13)
//   unvaccinated-stratum incidence: incV_u, incF_u, incN_u (families 14..16)
//   vaccinated-stratum incidence:   incV_v, incF_v, incN_v (families 17..19)
// Total = 20 * n_age = 80 states.


// ============================================================================

functions {
  // sum 14 carriage families per age to get pop_a from the state vector.
  // (families 14..19 are incidence accumulators not included)
  vector compute_pop_per_age(vector s, int n_age) {
    vector[n_age] pop_a = rep_vector(0.0, n_age);
    for (k in 0:13) {
      for (i in 1:n_age) {
        pop_a[i] += s[k * n_age + i];
      }
    }
    return pop_a;
  }

  // symmetrise daily contact matrix against supplied per-age pop and return annual contacts.
  matrix symmetrise_contacts(matrix C_raw_daily, vector pop_vec, int n_age) {
    matrix[n_age, n_age] Cs;
    for (i in 1:n_age) {
      for (j in 1:n_age) {
        Cs[i, j] = (C_raw_daily[i, j]
                  + C_raw_daily[j, i] * pop_vec[j] / pop_vec[i]) / 2.0;
      }
    }
    return Cs * 365.25;
  }

  // RHS of the within-year SIS ODE. Population per age (pop_j) and symmetrised contact matrix (c_sym) are passed in as constants
  // within-year carriage dynamics conserve pop_j.
  vector pneumo_rhs(real t, vector y,
                    int n_age,
                    vector q_age,
                    vector gamma_V, vector gamma_F, vector gamma_N,
                    real rr_V,  real rr_F,  real rr_N,
                    real rr_V7, real rr_F7, real rr_N7,
                    real delta_V, real delta_F,
                    matrix c_sym, vector pop_j) {
    int n  = n_age;
    int S0  = 0;       int V0  = 1 * n;   int F0  = 2 * n;   int N0  = 3 * n;
    int VF0 = 4 * n;   int NV0 = 5 * n;   int NF0 = 6 * n;
    int S70 = 7 * n;   int V70 = 8 * n;   int F70 = 9 * n;   int N70 = 10 * n;
    int VF70 = 11 * n; int NV70 = 12 * n; int NF70 = 13 * n;
    int IVu0 = 14 * n; int IFu0 = 15 * n; int INu0 = 16 * n;
    int IVv0 = 17 * n; int IFv0 = 18 * n; int INv0 = 19 * n;

    vector[n] S   = segment(y, S0  + 1, n);
    vector[n] V   = segment(y, V0  + 1, n);
    vector[n] Fc  = segment(y, F0  + 1, n);
    vector[n] Nc  = segment(y, N0  + 1, n);
    vector[n] VF  = segment(y, VF0 + 1, n);
    vector[n] NV  = segment(y, NV0 + 1, n);
    vector[n] NF  = segment(y, NF0 + 1, n);
    vector[n] S7  = segment(y, S70 + 1, n);
    vector[n] V7  = segment(y, V70 + 1, n);
    vector[n] F7  = segment(y, F70 + 1, n);
    vector[n] N7  = segment(y, N70 + 1, n);
    vector[n] VF7 = segment(y, VF70 + 1, n);
    vector[n] NV7 = segment(y, NV70 + 1, n);
    vector[n] NF7 = segment(y, NF70 + 1, n);

    vector[n] carrV = V  + VF  + NV  + V7  + VF7 + NV7;
    vector[n] carrF = Fc + VF  + NF  + F7  + VF7 + NF7;
    vector[n] carrN = Nc + NV  + NF  + N7  + NV7 + NF7;

    vector[n] lamV = q_age .* (c_sym * (carrV ./ pop_j));
    vector[n] lamF = q_age .* (c_sym * (carrF ./ pop_j));
    vector[n] lamN = q_age .* (c_sym * (carrN ./ pop_j));

    vector[n] lamV7 = (1.0 - delta_V) * lamV;
    vector[n] lamF7 = (1.0 - delta_F) * lamF;
    vector[n] lamN7 = lamN;

    vector[n] dS  = -S  .* (lamV + lamF + lamN) + gamma_V .* V  + gamma_F .* Fc + gamma_N .* Nc;
    vector[n] dV  =  S  .* lamV - V  .* (gamma_V + rr_V * (lamF + lamN)) + gamma_F .* VF + gamma_N .* NV;
    vector[n] dF  =  S  .* lamF - Fc .* (gamma_F + rr_F * (lamV + lamN)) + gamma_V .* VF + gamma_N .* NF;
    vector[n] dN  =  S  .* lamN - Nc .* (gamma_N + rr_N * (lamV + lamF)) + gamma_V .* NV + gamma_F .* NF;
    vector[n] dVF = rr_V * (V  .* lamF) + rr_F * (Fc .* lamV) - VF .* (gamma_V + gamma_F);
    vector[n] dNV = rr_V * (V  .* lamN) + rr_N * (Nc .* lamV) - NV .* (gamma_V + gamma_N);
    vector[n] dNF = rr_N * (Nc .* lamF) + rr_F * (Fc .* lamN) - NF .* (gamma_N + gamma_F);

    vector[n] dS7  = -S7  .* (lamV7 + lamF7 + lamN7) + gamma_V .* V7  + gamma_F .* F7 + gamma_N .* N7;
    vector[n] dV7  =  S7  .* lamV7 - V7  .* (gamma_V + rr_V7 * (lamF7 + lamN7)) + gamma_F .* VF7 + gamma_N .* NV7;
    vector[n] dF7  =  S7  .* lamF7 - F7  .* (gamma_F + rr_F7 * (lamV7 + lamN7)) + gamma_V .* VF7 + gamma_N .* NF7;
    vector[n] dN7  =  S7  .* lamN7 - N7  .* (gamma_N + rr_N7 * (lamV7 + lamF7)) + gamma_V .* NV7 + gamma_F .* NF7;
    vector[n] dVF7 = rr_V7 * (V7 .* lamF7) + rr_F7 * (F7 .* lamV7) - VF7 .* (gamma_V + gamma_F);
    vector[n] dNV7 = rr_V7 * (V7 .* lamN7) + rr_N7 * (N7 .* lamV7) - NV7 .* (gamma_V + gamma_N);
    vector[n] dNF7 = rr_N7 * (N7 .* lamF7) + rr_F7 * (F7 .* lamN7) - NF7 .* (gamma_N + gamma_F);

    vector[n] d_incV_u = S  .* lamV  + rr_F  * (Fc .* lamV)  + rr_N  * (Nc .* lamV);
    vector[n] d_incF_u = S  .* lamF  + rr_V  * (V  .* lamF)  + rr_N  * (Nc .* lamF);
    vector[n] d_incN_u = S  .* lamN  + rr_V  * (V  .* lamN)  + rr_F  * (Fc .* lamN);

    vector[n] d_incV_v = S7 .* lamV7 + rr_F7 * (F7 .* lamV7) + rr_N7 * (N7 .* lamV7);
    vector[n] d_incF_v = S7 .* lamF7 + rr_V7 * (V7 .* lamF7) + rr_N7 * (N7 .* lamF7);
    vector[n] d_incN_v = S7 .* lamN7 + rr_V7 * (V7 .* lamN7) + rr_F7 * (F7 .* lamN7);

    int n_states = 20 * n;
    vector[n_states] dy = rep_vector(0.0, n_states);
    for (i in 1:n) {
      dy[S0  + i] = dS[i];   dy[V0  + i] = dV[i];   dy[F0  + i] = dF[i];
      dy[N0  + i] = dN[i];   dy[VF0 + i] = dVF[i];  dy[NV0 + i] = dNV[i];
      dy[NF0 + i] = dNF[i];
      dy[S70 + i] = dS7[i];  dy[V70 + i] = dV7[i];  dy[F70 + i] = dF7[i];
      dy[N70 + i] = dN7[i];  dy[VF70 + i] = dVF7[i]; dy[NV70 + i] = dNV7[i];
      dy[NF70 + i] = dNF7[i];
      dy[IVu0 + i] = d_incV_u[i]; dy[IFu0 + i] = d_incF_u[i]; dy[INu0 + i] = d_incN_u[i];
      dy[IVv0 + i] = d_incV_v[i]; dy[IFv0 + i] = d_incF_v[i]; dy[INv0 + i] = d_incN_v[i];
    }
    return dy;
  }

  // year-boundary update for pre-PCV burn-in: waning + aging + birth injection equal to target_pop[1] + rescale-to-target_pop + reset incidence
  // rescale keeps total pop at init_pop_1999 throughout burn-in so carriage state can settle to its 1999-level equilibrium.
  vector age_rescale_burnin(vector y, int n_age, vector aging_rate,
                            real vacc_cov, real omega_V, real omega_F,
                            vector target_pop) {
    int n = n_age;
    int n_states = 20 * n;
    vector[n_states] s = y;
    real wfV = 1.0 - exp(-omega_V);
    real wfF = 1.0 - exp(-omega_F);

    // waning of PCV7 protection
    for (i in 1:n) {
      int S_  = i;            int V_  = 1*n + i;  int F_  = 2*n + i;
      int N_  = 3*n + i;      int VF_ = 4*n + i;  int NV_ = 5*n + i;
      int NF_ = 6*n + i;
      int S7_ = 7*n + i;      int V7_ = 8*n + i;  int F7_ = 9*n + i;
      int N7_ = 10*n + i;     int VF7_= 11*n + i; int NV7_= 12*n + i;
      int NF7_= 13*n + i;
      real m_S  = s[S7_]  * wfV; s[S7_]  -= m_S;  s[S_]  += m_S;
      real m_V  = s[V7_]  * wfV; s[V7_]  -= m_V;  s[V_]  += m_V;
      real m_N  = s[N7_]  * wfV; s[N7_]  -= m_N;  s[N_]  += m_N;
      real m_NV = s[NV7_] * wfV; s[NV7_] -= m_NV; s[NV_] += m_NV;
      real n_F  = s[F7_]  * wfF; s[F7_]  -= n_F;  s[F_]  += n_F;
      real n_VF = s[VF7_] * wfF; s[VF7_] -= n_VF; s[VF_] += n_VF;
      real n_NF = s[NF7_] * wfF; s[NF7_] -= n_NF; s[NF_] += n_NF;
    }

    // aging (carriage families 0..13 only)
    for (k in 0:13) {
      int base = k * n;
      vector[n] orig;
      for (i in 1:n) orig[i] = s[base + i];
      for (i in 1:n) {
        real new_i = orig[i] * (1.0 - aging_rate[i]);
        if (i > 1) new_i += orig[i - 1] * aging_rate[i - 1];
        s[base + i] = new_i;
      }
    }

    // births equal to target_pop[1] (stationary burn-in)
    {
      real to_unvacc = (1.0 - vacc_cov) * target_pop[1];
      real to_vacc   =        vacc_cov  * target_pop[1];
      s[1]           = s[1]         + to_unvacc;
      s[7 * n + 1]   = s[7 * n + 1] + to_vacc;
    }

    // rescale each age to target_pop
    for (i in 1:n) {
      real total = 0.0;
      for (k in 0:13) total += s[k * n + i];
      if (total > 0) {
        real scale = target_pop[i] / total;
        for (k in 0:13) s[k * n + i] *= scale;
      }
    }

    // reset incidence accumulators
    for (k in 14:19) for (i in 1:n) s[k * n + i] = 0.0;
    return s;
  }

  // year-boundary update for the fit loop: waning, then a 12-substep Euler integration of the demographic ODE for one year, then reset incidence.
  // no rescale, population evolves freely under explicit birth/death/aging.
  //
  // demographic ODE per substep (dt = 1/12 yr), applied UNIFORMLY across all 14 carriage families per age (carriage state is preserved):
  //   For each age a:
  //     out_rate_a = aging_rate[a] + mu_age[a]
  //     N_a(t + dt) = N_a(t) * (1 - out_rate_a * dt) + N_{a-1}(t) * aging_rate[a-1] * dt    (for a > 1)
  //   births at the substep = birth_rate_y * Total(t) * dt
  //   births split between S<1y and S_7<1y by vacc_cov.
  // empties into 1-4y in one year, refilled only by births.
  vector age_demog_step(vector y, int n_age, vector aging_rate, vector mu_age,
                        real vacc_cov, real omega_V, real omega_F,
                        real birth_rate_y) {
    int n = n_age;
    int n_states = 20 * n;
    vector[n_states] s = y;
    real wfV = 1.0 - exp(-omega_V);
    real wfF = 1.0 - exp(-omega_F);
    int n_substeps = 12;
    real dt = 1.0 / n_substeps;

    // waning of PCV7 protection (instantaneous, at year boundary)
    for (i in 1:n) {
      int S_  = i;            int V_  = 1*n + i;  int F_  = 2*n + i;
      int N_  = 3*n + i;      int VF_ = 4*n + i;  int NV_ = 5*n + i;
      int NF_ = 6*n + i;
      int S7_ = 7*n + i;      int V7_ = 8*n + i;  int F7_ = 9*n + i;
      int N7_ = 10*n + i;     int VF7_= 11*n + i; int NV7_= 12*n + i;
      int NF7_= 13*n + i;
      real m_S  = s[S7_]  * wfV; s[S7_]  -= m_S;  s[S_]  += m_S;
      real m_V  = s[V7_]  * wfV; s[V7_]  -= m_V;  s[V_]  += m_V;
      real m_N  = s[N7_]  * wfV; s[N7_]  -= m_N;  s[N_]  += m_N;
      real m_NV = s[NV7_] * wfV; s[NV7_] -= m_NV; s[NV_] += m_NV;
      real n_F  = s[F7_]  * wfF; s[F7_]  -= n_F;  s[F_]  += n_F;
      real n_VF = s[VF7_] * wfF; s[VF7_] -= n_VF; s[VF_] += n_VF;
      real n_NF = s[NF7_] * wfF; s[NF7_] -= n_NF; s[NF_] += n_NF;
    }

    // demographic step: mortality + aging (uniform across 14 carriage families per age) + births (into S<1y / S_7<1y only) over one year, integrated in n_substeps Euler steps.
    for (sub in 1:n_substeps) {
      // per-age total pop = sum over 14 carriage families
      vector[n] pop_a = compute_pop_per_age(s, n);
      real total_pop = sum(pop_a);
      real births    = birth_rate_y * total_pop;

      // mortality + aging applied to every carriage family
      for (k in 0:13) {
        int base = k * n;
        vector[n] orig;
        for (i in 1:n) orig[i] = s[base + i];
        for (i in 1:n) {
          real out_rate = aging_rate[i] + mu_age[i];
          real new_i    = orig[i] * (1.0 - out_rate * dt);
          if (i > 1) new_i += orig[i - 1] * aging_rate[i - 1] * dt;
          s[base + i] = new_i;
        }
      }

      // births enter S<1y / S_7<1y only
      s[1]         += (1.0 - vacc_cov) * births * dt;
      s[7 * n + 1] +=        vacc_cov  * births * dt;
    }

    // reset incidence accumulators (families 14..19)
    for (k in 14:19) for (i in 1:n) s[k * n + i] = 0.0;
    return s;
  }
}

// ============================================================================

data {
  int<lower=1> n_age; // 4
  int<lower=1> n_year; // 1999..2009 = 11
  int<lower=1, upper=n_year> n_fit;   // first n_fit years enter the Poisson likelihood; ALL_YEARS[(n_fit+1):n_year] are out-of-sample projections
  int<lower=1> n_burn; // burn-in years inside Stan
  
  vector[n_age]              gamma_V;
  vector[n_age]              gamma_F;
  vector[n_age]              gamma_N;
  vector[n_age]              aging_rate;

  // unsymmetrised daily contact matrix, symmetrised inside Stan against current model populations.
  matrix[n_age, n_age]       contact_raw;

  // initial 1999 population (calibrated by 06_demographic_calibration.R) used as burn-in rescale target, seed state's totals, and year-1 (1999) FOI normalisation.
  vector<lower=0>[n_age]     init_pop_1999;

  // explicit demographic inputs. birth_rate_year[t] is crude per-capita birth rate (per year) for year t of fit window. 
  // mu_age[a] is age-specific net mortality rate per year (if negative represents mortality - net immigration).
  vector<lower=0>[n_year]    birth_rate_year;
  vector[n_age]              mu_age;

  // seed state from R lightweight burn-in (length 20 * n_age).
  vector[20 * n_age]         seed_state;

  // independent log-uniform bounds (per-parameter) on the log-scale uncertainty-propagation parameters. 
  vector[n_age]              log_q_lower;
  vector[n_age]              log_q_upper;
  real                       log_rr_V_lb; real log_rr_V_ub;
  real                       log_rr_F_lb; real log_rr_F_ub;
  real                       log_rr_N_lb; real log_rr_N_ub;

  // prior 95 % CrI bounds for delta and omega (qbeta/qgamma in R).
  real<lower=0, upper=1>     delta_prior_lo;
  real<lower=0, upper=1>     delta_prior_hi;
  real<lower=0>              omega_prior_lo;
  real<lower=0>              omega_prior_hi;

  array[n_year, n_age, 3] int<lower=0> obs_ipd;

  // per-boundary vaccination-at-birth coverage; current schedule
  // (0.43, 0.98, 0.98, 0.98, 0.95, ..., 0.95).
  vector<lower=0, upper=1>[n_year] vacc_cov_year;
}

// ============================================================================

transformed data {
  int n_states = 20 * n_age;
  real t0      = 0.0;
  array[1] real ts = {1.0};

  real ode_rel_tol     = 1e-4;
  real ode_abs_tol     = 1e-4;
  int  ode_max_num_steps = 10000;

  // log-scale bounds for estimated parameters tied to prior 95 % CrI.
  real log_delta_lb = log(delta_prior_lo);
  real log_delta_ub = log(delta_prior_hi);
  real log_omega_lb = log(omega_prior_lo);
  real log_omega_ub = log(omega_prior_hi);

  // c_sym at 1999 levels, used during the in-Stan burn-in (which keeps pop at init_pop_1999).
  matrix[n_age, n_age] c_sym_1999 = symmetrise_contacts(contact_raw,
                                                        init_pop_1999, n_age);
}

// ============================================================================

parameters {
  // log-scale estimated parameters, constraints tied to prior 95 % CrI.
  real<lower=log_delta_lb, upper=log_delta_ub>  log_delta_V;
  real<lower=log_delta_lb, upper=log_delta_ub>  log_delta_F;
  real<lower=log_omega_lb, upper=log_omega_ub>  log_omega_V;
  real<lower=log_omega_lb, upper=log_omega_ub>  log_omega_F;

  // log-scale uncertainty-propagation parameters with independent log-uniform priors (bounds = marginal quantiles of the carriage-fit posterior).
  // Stan's default uniform on unconstrained scale, with these <lower=, upper=> constraints and no explicit prior in the model block, gives log-uniform on natural scale.
  vector<lower=log_q_lower, upper=log_q_upper>[n_age] log_q_age;
  real<lower=log_rr_V_lb, upper=log_rr_V_ub> log_rr_V;
  real<lower=log_rr_F_lb, upper=log_rr_F_ub> log_rr_F;
  real<lower=log_rr_N_lb, upper=log_rr_N_ub> log_rr_N;
}

// ============================================================================

transformed parameters {
  // natural-scale estimated parameters
  real delta_V = exp(log_delta_V);
  real delta_F = exp(log_delta_F);
  real omega_V = exp(log_omega_V);
  real omega_F = exp(log_omega_F);

  // natural-scale uncertainty-propagation parameters, log_* are sampled directly with bounded constraints, recovered by exp().
  vector[n_age] q_age = exp(log_q_age);
  real rr_V = exp(log_rr_V);
  real rr_F = exp(log_rr_F);
  real rr_N = exp(log_rr_N);
  real rr_V7 = rr_V;
  real rr_F7 = rr_F;
  real rr_N7 = rr_N;

  // per-year arrays
  array[n_year, n_age, 3] real pred_ipd;
  array[n_year, n_age, 3] real inc_year;
  array[n_year, n_age, 3] real inc_year_unvacc;
  array[n_year, n_age, 3] real inc_year_vacc;
  array[n_year] vector[n_states] year_end_state;
  array[n_year] vector[n_states] year_start_state;
  array[n_year] vector[n_age]    pop_year_start; // model-derived pop per age
  matrix[n_age, 3] ccr_1999;

  vector[n_states] cur = seed_state;

  // pre-PCV burn-in: stationary at init_pop_1999, no demographics (carriage state equilibrates; pop totals pinned by rescale)
  for (b in 1:n_burn) {
    array[1] vector[n_states] sol_b = ode_bdf_tol(pneumo_rhs, cur, t0, ts,
                                               ode_rel_tol, ode_abs_tol,
                                               ode_max_num_steps,
                                               n_age, q_age,
                                               gamma_V, gamma_F, gamma_N,
                                               rr_V, rr_F, rr_N,
                                               rr_V7, rr_F7, rr_N7,
                                               delta_V, delta_F,
                                               c_sym_1999, init_pop_1999);
    cur = age_rescale_burnin(sol_b[1], n_age, aging_rate,
                              0.0, omega_V, omega_F, init_pop_1999);
  }

  // main fit loop, years 1999..2009 with explicit demographics
  for (t in 1:n_year) {
    year_start_state[t] = cur;

    // compute model populations at year start from compartment sums and symmetrise contact matrix against them
    // FOI normalisation and contact intensities are internally consistent with carriage state.
    vector[n_age] pop_now = compute_pop_per_age(cur, n_age);
    pop_year_start[t]     = pop_now;
    matrix[n_age, n_age] c_sym_now = symmetrise_contacts(contact_raw,
                                                         pop_now, n_age);

    array[1] vector[n_states] sol = ode_bdf_tol(pneumo_rhs, cur, t0, ts,
                                             ode_rel_tol, ode_abs_tol,
                                             ode_max_num_steps,
                                             n_age, q_age,
                                             gamma_V, gamma_F, gamma_N,
                                             rr_V, rr_F, rr_N,
                                             rr_V7, rr_F7, rr_N7,
                                             delta_V, delta_F,
                                             c_sym_now, pop_now);
    vector[n_states] endv = sol[1];

    for (a in 1:n_age) {
      inc_year_unvacc[t, a, 1] = endv[14 * n_age + a];
      inc_year_unvacc[t, a, 2] = endv[15 * n_age + a];
      inc_year_unvacc[t, a, 3] = endv[16 * n_age + a];
      inc_year_vacc[t, a, 1]   = endv[17 * n_age + a];
      inc_year_vacc[t, a, 2]   = endv[18 * n_age + a];
      inc_year_vacc[t, a, 3]   = endv[19 * n_age + a];
      inc_year[t, a, 1] = inc_year_unvacc[t, a, 1] + inc_year_vacc[t, a, 1];
      inc_year[t, a, 2] = inc_year_unvacc[t, a, 2] + inc_year_vacc[t, a, 2];
      inc_year[t, a, 3] = inc_year_unvacc[t, a, 3] + inc_year_vacc[t, a, 3];
    }

    if (t == 1) {
      for (a in 1:n_age) {
        for (k in 1:3) {
          real inc_clip = fmax(inc_year[t, a, k], 1e-6);
          ccr_1999[a, k] = obs_ipd[t, a, k] / inc_clip;
        }
      }
    }

    for (a in 1:n_age) {
      pred_ipd[t, a, 1] = (inc_year_unvacc[t, a, 1] + inc_year_vacc[t, a, 1]) * ccr_1999[a, 1];
      pred_ipd[t, a, 2] = (inc_year_unvacc[t, a, 2] + inc_year_vacc[t, a, 2]) * ccr_1999[a, 2];
      pred_ipd[t, a, 3] = (inc_year_unvacc[t, a, 3] + inc_year_vacc[t, a, 3]) * ccr_1999[a, 3];
    }

    year_end_state[t] = endv;

    // explicit demographic year-boundary update.
    real vacc_now = vacc_cov_year[t];
    cur = age_demog_step(endv, n_age, aging_rate, mu_age,
                          vacc_now, omega_V, omega_F,
                          birth_rate_year[t]);
  }
}

// ============================================================================

model {
  delta_V ~ beta(13, 7);
  delta_F ~ beta(13, 7);
  omega_V ~ gamma(4, 11);
  omega_F ~ gamma(4, 11);
  target += log_delta_V + log_delta_F + log_omega_V + log_omega_F;

  // q_age and rr_V/F/N have independent log-uniform priors implied by <lower=, upper=> constraints 
  // only the first n_fit years enter the likelihood. Years (n_fit+1):n_year still have their dynamics simulated in transformed parameters
  // pred_ipd / pred_ipd_rep are generated for them (out-of-sample projections), but the observed IPD counts for those years do not contribute to the posterior.
  for (t in 1:n_fit) {
    for (a in 1:n_age) {
      for (k in 1:3) {
        target += poisson_lpmf(obs_ipd[t, a, k] | fmax(pred_ipd[t, a, k], 0.0) + 1e-6);
      }
    }
  }
}

// ============================================================================

generated quantities {
  array[n_year, n_age, 3] real foi_year;
  array[n_year, n_age, 3] real carr_year_start;
  array[n_year, n_age, 3] real prev_year_start;
  array[n_year, n_age, 3] real carr_year_end;
  array[n_year, n_age, 3] real prev_year_end;
  array[n_year, n_age]    real model_pop_year_start;  // total pop per age, year start
  array[n_year, n_age]    real model_pop_year_end;    // total pop per age, year end
  array[n_year, n_age, 3] int  pred_ipd_rep;
  array[n_year, n_age, 3] real log_lik;
  real duration_V = 1.0 / omega_V;
  real duration_F = 1.0 / omega_F;

  // FOI uses model populations and symmetrised-against-model contacts.
  for (t in 1:n_year) {
    vector[n_age] carrV;
    vector[n_age] carrF;
    vector[n_age] carrN;
    vector[n_states] yss = year_start_state[t];
    for (i in 1:n_age) {
      carrV[i] = yss[1*n_age + i] + yss[4*n_age + i] + yss[5*n_age + i]
                + yss[8*n_age + i] + yss[11*n_age + i] + yss[12*n_age + i];
      carrF[i] = yss[2*n_age + i] + yss[4*n_age + i] + yss[6*n_age + i]
                + yss[9*n_age + i] + yss[11*n_age + i] + yss[13*n_age + i];
      carrN[i] = yss[3*n_age + i] + yss[5*n_age + i] + yss[6*n_age + i]
                + yss[10*n_age + i] + yss[12*n_age + i] + yss[13*n_age + i];
    }
    vector[n_age] pop_s = pop_year_start[t];
    matrix[n_age, n_age] c_sym_s = symmetrise_contacts(contact_raw, pop_s, n_age);
    vector[n_age] lamV = q_age .* (c_sym_s * (carrV ./ pop_s));
    vector[n_age] lamF = q_age .* (c_sym_s * (carrF ./ pop_s));
    vector[n_age] lamN = q_age .* (c_sym_s * (carrN ./ pop_s));

    vector[n_age] carrV_e;
    vector[n_age] carrF_e;
    vector[n_age] carrN_e;
    vector[n_states] yse = year_end_state[t];
    for (i in 1:n_age) {
      carrV_e[i] = yse[1*n_age + i] + yse[4*n_age + i] + yse[5*n_age + i]
                  + yse[8*n_age + i] + yse[11*n_age + i] + yse[12*n_age + i];
      carrF_e[i] = yse[2*n_age + i] + yse[4*n_age + i] + yse[6*n_age + i]
                  + yse[9*n_age + i] + yse[11*n_age + i] + yse[13*n_age + i];
      carrN_e[i] = yse[3*n_age + i] + yse[5*n_age + i] + yse[6*n_age + i]
                  + yse[10*n_age + i] + yse[12*n_age + i] + yse[13*n_age + i];
    }
    vector[n_age] pop_e = compute_pop_per_age(yse, n_age);

    for (a in 1:n_age) {
      foi_year[t, a, 1] = lamV[a];
      foi_year[t, a, 2] = lamF[a];
      foi_year[t, a, 3] = lamN[a];
      carr_year_start[t, a, 1] = carrV[a];
      carr_year_start[t, a, 2] = carrF[a];
      carr_year_start[t, a, 3] = carrN[a];
      prev_year_start[t, a, 1] = carrV[a] / pop_s[a];
      prev_year_start[t, a, 2] = carrF[a] / pop_s[a];
      prev_year_start[t, a, 3] = carrN[a] / pop_s[a];
      carr_year_end[t, a, 1] = carrV_e[a];
      carr_year_end[t, a, 2] = carrF_e[a];
      carr_year_end[t, a, 3] = carrN_e[a];
      prev_year_end[t, a, 1] = carrV_e[a] / pop_e[a];
      prev_year_end[t, a, 2] = carrF_e[a] / pop_e[a];
      prev_year_end[t, a, 3] = carrN_e[a] / pop_e[a];
      model_pop_year_start[t, a] = pop_s[a];
      model_pop_year_end[t, a]   = pop_e[a];
      for (k in 1:3) {
        real lambda_k = fmax(pred_ipd[t, a, k], 0.0) + 1e-6;
        pred_ipd_rep[t, a, k] = poisson_rng(lambda_k);
        log_lik[t, a, k]      = poisson_lpmf(obs_ipd[t, a, k] | lambda_k);
      }
    }
  }
}
