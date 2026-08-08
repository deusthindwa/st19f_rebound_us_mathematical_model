// by Deus Thindwa, Dan Weinberger, Ginny Pitzer
// age-structured mathematical model for pneumococcal transmission
// 05/23/2026

// =============================================================================

// PCV7 -> PCV13 model with JOINT-POSTERIOR (multivariate-normal) propagation.
// This code serves all three PCV13 scenarios.
//
//   - Initial 1999 population from `init_1999_fitted` (demographic_model.rds).
//   - year-boundary update applies explicit births and per-age mortality as pre-pcv13
//   - propagated parameters are drawn from two MULTIVARIATE-NORMAL priors on the log scale,
//     with mean vector and Cholesky-factored covariance matrix computed in
//     04_fit_pcv13_common.R from the JOINT posteriors of the source fits:
//        carriage_fit -> log_carr = [log_q_age(1..n_age), log_rr_V, log_rr_F, log_rr_N_pre]
//        ipd_fit      -> log_ipd  = [log_delta_V_7, log_delta_F_7, log_omega_V_7, log_omega_F_7]
//     Using multi_normal_cholesky preserves ALL pairwise correlations from the
//     source posteriors (under a Gaussian assumption on log scale).
//   - delta_V_13 = delta_V_7 and omega_V_13 = omega_V_7
//   - rr_N is period-specific: rr_N_pre drives 1999-2009 dynamics
//   - rr_N_post drives 2010-2019 dynamics. rr_N_post is estimated in scenario 3, prior Beta(7, 2), for scenarios 1 and 2 it is sampled the same way as rr_N_pre.
//   - scenario-specific informative priors (Beta / Gamma) with Jacobian:
//        scenario 1: delta_F_13 ~ Beta(3, 9)
//        scenario 2: omega_F_13 ~ Gamma(shape = 3, rate = 3)
//        scenario 3: rr_N_post ~ Beta(7, 2)
//   - Likelihood summed only over PCV13 fit window 2010..2019.
//
// state vector layout per age (24 families = 21 carriage + 3 incidence):
//   Unvaccinated (0..6):    S, V, F, N, VF, NV, NF
//   PCV7-vaccinated (7..13): S_7, V_7, F_7, N_7, VF_7, NV_7, NF_7
//   PCV13-vaccinated (14..20): S_13, V_13, F_13, N_13, VF_13, NV_13, NF_13
//   Total incidence (21..23): incV, incF, incN
// Total = 24 * n_age = 96 states for n_age = 4.

// =============================================================================

functions {
  // Within-year SIS ODE (3 strata). pop_j is the model population per age
  // at the start of the year, used in the FOI denominator and contact-matrix symmetrisation
  vector pneumo_rhs(real t, vector y,
                    int n_age,
                    vector q_age,
                    vector gamma_V, vector gamma_F, vector gamma_N,
                    real rr_V, real rr_F, real rr_N,
                    real delta_V_7, real delta_F_7,
                    real delta_V_13, real delta_F_13,
                    matrix c_sym, vector pop_j) {
    int n = n_age;
    int S0   = 0;        int V0   = 1*n;   int F0   = 2*n;   int N0   = 3*n;
    int VF0  = 4*n;      int NV0  = 5*n;   int NF0  = 6*n;
    int S70  = 7*n;      int V70  = 8*n;   int F70  = 9*n;   int N70  = 10*n;
    int VF70 = 11*n;     int NV70 = 12*n;  int NF70 = 13*n;
    int S130 = 14*n;     int V130 = 15*n;  int F130 = 16*n;  int N130 = 17*n;
    int VF130= 18*n;     int NV130= 19*n;  int NF130= 20*n;
    int IV0  = 21*n;     int IF0  = 22*n;  int IN0  = 23*n;

    vector[n] S    = segment(y, S0   + 1, n);
    vector[n] V    = segment(y, V0   + 1, n);
    vector[n] Fc   = segment(y, F0   + 1, n);
    vector[n] Nc   = segment(y, N0   + 1, n);
    vector[n] VF   = segment(y, VF0  + 1, n);
    vector[n] NV   = segment(y, NV0  + 1, n);
    vector[n] NF   = segment(y, NF0  + 1, n);
    vector[n] S7   = segment(y, S70  + 1, n);
    vector[n] V7   = segment(y, V70  + 1, n);
    vector[n] F7   = segment(y, F70  + 1, n);
    vector[n] N7   = segment(y, N70  + 1, n);
    vector[n] VF7  = segment(y, VF70 + 1, n);
    vector[n] NV7  = segment(y, NV70 + 1, n);
    vector[n] NF7  = segment(y, NF70 + 1, n);
    vector[n] S13  = segment(y, S130 + 1, n);
    vector[n] V13  = segment(y, V130 + 1, n);
    vector[n] F13  = segment(y, F130 + 1, n);
    vector[n] N13  = segment(y, N130 + 1, n);
    vector[n] VF13 = segment(y, VF130+ 1, n);
    vector[n] NV13 = segment(y, NV130+ 1, n);
    vector[n] NF13 = segment(y, NF130+ 1, n);

    vector[n] carrV = V  + VF  + NV  + V7  + VF7  + NV7  + V13  + VF13  + NV13;
    vector[n] carrF = Fc + VF  + NF  + F7  + VF7  + NF7  + F13  + VF13  + NF13;
    vector[n] carrN = Nc + NV  + NF  + N7  + NV7  + NF7  + N13  + NV13  + NF13;

    vector[n] lamV = q_age .* (c_sym * (carrV ./ pop_j));
    vector[n] lamF = q_age .* (c_sym * (carrF ./ pop_j));
    vector[n] lamN = q_age .* (c_sym * (carrN ./ pop_j));

    vector[n] lamV_7  = (1.0 - delta_V_7)  * lamV;
    vector[n] lamF_7  = (1.0 - delta_F_7)  * lamF;
    vector[n] lamN_7  = lamN;
    vector[n] lamV_13 = (1.0 - delta_V_13) * lamV;
    vector[n] lamF_13 = (1.0 - delta_F_13) * lamF;
    vector[n] lamN_13 = lamN;

    vector[n] dS  = -S  .* (lamV + lamF + lamN) + gamma_V .* V  + gamma_F .* Fc + gamma_N .* Nc;
    vector[n] dV  =  S  .* lamV - V  .* (gamma_V + rr_V * (lamF + lamN)) + gamma_F .* VF + gamma_N .* NV;
    vector[n] dF  =  S  .* lamF - Fc .* (gamma_F + rr_F * (lamV + lamN)) + gamma_V .* VF + gamma_N .* NF;
    vector[n] dN  =  S  .* lamN - Nc .* (gamma_N + rr_N * (lamV + lamF)) + gamma_V .* NV + gamma_F .* NF;
    vector[n] dVF = rr_V * (V  .* lamF) + rr_F * (Fc .* lamV) - VF .* (gamma_V + gamma_F);
    vector[n] dNV = rr_V * (V  .* lamN) + rr_N * (Nc .* lamV) - NV .* (gamma_V + gamma_N);
    vector[n] dNF = rr_N * (Nc .* lamF) + rr_F * (Fc .* lamN) - NF .* (gamma_N + gamma_F);

    vector[n] dS7  = -S7  .* (lamV_7 + lamF_7 + lamN_7) + gamma_V .* V7  + gamma_F .* F7 + gamma_N .* N7;
    vector[n] dV7  =  S7  .* lamV_7 - V7  .* (gamma_V + rr_V * (lamF_7 + lamN_7)) + gamma_F .* VF7 + gamma_N .* NV7;
    vector[n] dF7  =  S7  .* lamF_7 - F7  .* (gamma_F + rr_F * (lamV_7 + lamN_7)) + gamma_V .* VF7 + gamma_N .* NF7;
    vector[n] dN7  =  S7  .* lamN_7 - N7  .* (gamma_N + rr_N * (lamV_7 + lamF_7)) + gamma_V .* NV7 + gamma_F .* NF7;
    vector[n] dVF7 = rr_V * (V7 .* lamF_7) + rr_F * (F7 .* lamV_7) - VF7 .* (gamma_V + gamma_F);
    vector[n] dNV7 = rr_V * (V7 .* lamN_7) + rr_N * (N7 .* lamV_7) - NV7 .* (gamma_V + gamma_N);
    vector[n] dNF7 = rr_N * (N7 .* lamF_7) + rr_F * (F7 .* lamN_7) - NF7 .* (gamma_N + gamma_F);

    vector[n] dS13  = -S13  .* (lamV_13 + lamF_13 + lamN_13) + gamma_V .* V13  + gamma_F .* F13 + gamma_N .* N13;
    vector[n] dV13  =  S13  .* lamV_13 - V13  .* (gamma_V + rr_V * (lamF_13 + lamN_13)) + gamma_F .* VF13 + gamma_N .* NV13;
    vector[n] dF13  =  S13  .* lamF_13 - F13  .* (gamma_F + rr_F * (lamV_13 + lamN_13)) + gamma_V .* VF13 + gamma_N .* NF13;
    vector[n] dN13  =  S13  .* lamN_13 - N13  .* (gamma_N + rr_N * (lamV_13 + lamF_13)) + gamma_V .* NV13 + gamma_F .* NF13;
    vector[n] dVF13 = rr_V * (V13 .* lamF_13) + rr_F * (F13 .* lamV_13) - VF13 .* (gamma_V + gamma_F);
    vector[n] dNV13 = rr_V * (V13 .* lamN_13) + rr_N * (N13 .* lamV_13) - NV13 .* (gamma_V + gamma_N);
    vector[n] dNF13 = rr_N * (N13 .* lamF_13) + rr_F * (F13 .* lamN_13) - NF13 .* (gamma_N + gamma_F);

    vector[n] d_incV =  S   .* lamV    + rr_F * (Fc  .* lamV)    + rr_N * (Nc  .* lamV)
                     +  S7  .* lamV_7  + rr_F * (F7  .* lamV_7)  + rr_N * (N7  .* lamV_7)
                     +  S13 .* lamV_13 + rr_F * (F13 .* lamV_13) + rr_N * (N13 .* lamV_13);
    vector[n] d_incF =  S   .* lamF    + rr_V * (V   .* lamF)    + rr_N * (Nc  .* lamF)
                     +  S7  .* lamF_7  + rr_V * (V7  .* lamF_7)  + rr_N * (N7  .* lamF_7)
                     +  S13 .* lamF_13 + rr_V * (V13 .* lamF_13) + rr_N * (N13 .* lamF_13);
    vector[n] d_incN =  S   .* lamN    + rr_V * (V   .* lamN)    + rr_F * (Fc  .* lamN)
                     +  S7  .* lamN_7  + rr_V * (V7  .* lamN_7)  + rr_F * (F7  .* lamN_7)
                     +  S13 .* lamN_13 + rr_V * (V13 .* lamN_13) + rr_F * (F13 .* lamN_13);

    int n_states = 24 * n;
    vector[n_states] dy = rep_vector(0.0, n_states);
    for (i in 1:n) {
      dy[S0   + i] = dS[i];   dy[V0   + i] = dV[i];   dy[F0   + i] = dF[i];
      dy[N0   + i] = dN[i];   dy[VF0  + i] = dVF[i];  dy[NV0  + i] = dNV[i];
      dy[NF0  + i] = dNF[i];
      dy[S70  + i] = dS7[i];  dy[V70  + i] = dV7[i];  dy[F70  + i] = dF7[i];
      dy[N70  + i] = dN7[i];  dy[VF70 + i] = dVF7[i]; dy[NV70 + i] = dNV7[i];
      dy[NF70 + i] = dNF7[i];
      dy[S130 + i] = dS13[i]; dy[V130 + i] = dV13[i]; dy[F130 + i] = dF13[i];
      dy[N130 + i] = dN13[i]; dy[VF130+ i] = dVF13[i];dy[NV130+ i] = dNV13[i];
      dy[NF130+ i] = dNF13[i];
      dy[IV0  + i] = d_incV[i]; dy[IF0  + i] = d_incF[i]; dy[IN0  + i] = d_incN[i];
    }
    return dy;
  }

  // Year-boundary update with
  //   1. waning of PCV7 and PCV13 protection.
  //   2. aging across age groups + per-age net mortality applied together via single Euler step over 1 year as pre-pcv13
  //   3. Vaccination at birth: total newborns = birth_rate * total_pop, split between unvacc, PCV7, PCV13 by vacc_cov_pcv7 / vacc_cov_pcv13.
  //   4. Reset incidence accumulators.
  vector age_vaccinate(vector y, int n_age, vector aging_rate, vector mu_age,
                       real vacc_cov_pcv7, real vacc_cov_pcv13,
                       real omega_V_7, real omega_F_7,
                       real omega_V_13, real omega_F_13,
                       real birth_rate) {
    int n = n_age;
    int n_states = 24 * n;
    vector[n_states] s = y;
    real wfV_7  = 1.0 - exp(-omega_V_7);
    real wfF_7  = 1.0 - exp(-omega_F_7);
    real wfV_13 = 1.0 - exp(-omega_V_13);
    real wfF_13 = 1.0 - exp(-omega_F_13);

    // waning (per age, PCV7 then PCV13)
    for (i in 1:n) {
      int S_   = i;            int V_   = 1*n + i;  int F_   = 2*n + i;
      int N_   = 3*n + i;      int VF_  = 4*n + i;  int NV_  = 5*n + i;
      int NF_  = 6*n + i;
      int S7_  = 7*n + i;      int V7_  = 8*n + i;  int F7_  = 9*n + i;
      int N7_  = 10*n + i;     int VF7_ = 11*n + i; int NV7_ = 12*n + i;
      int NF7_ = 13*n + i;
      int S13_ = 14*n + i;     int V13_ = 15*n + i; int F13_ = 16*n + i;
      int N13_ = 17*n + i;     int VF13_= 18*n + i; int NV13_= 19*n + i;
      int NF13_= 20*n + i;
      // PCV7 V-protected pool to unvacc
      real m1 = s[S7_]  * wfV_7;  s[S7_]  -= m1; s[S_]  += m1;
      real m2 = s[V7_]  * wfV_7;  s[V7_]  -= m2; s[V_]  += m2;
      real m3 = s[N7_]  * wfV_7;  s[N7_]  -= m3; s[N_]  += m3;
      real m4 = s[NV7_] * wfV_7;  s[NV7_] -= m4; s[NV_] += m4;
      // PCV7 F-protected pool to unvacc
      real m5 = s[F7_]  * wfF_7;  s[F7_]  -= m5; s[F_]  += m5;
      real m6 = s[VF7_] * wfF_7;  s[VF7_] -= m6; s[VF_] += m6;
      real m7 = s[NF7_] * wfF_7;  s[NF7_] -= m7; s[NF_] += m7;
      // PCV13 V-protected pool to unvacc
      real m8  = s[S13_]  * wfV_13; s[S13_]  -= m8;  s[S_]  += m8;
      real m9  = s[V13_]  * wfV_13; s[V13_]  -= m9;  s[V_]  += m9;
      real m10 = s[N13_]  * wfV_13; s[N13_]  -= m10; s[N_]  += m10;
      real m11 = s[NV13_] * wfV_13; s[NV13_] -= m11; s[NV_] += m11;
      // PCV13 F-protected pool to unvacc
      real m12 = s[F13_]  * wfF_13; s[F13_]  -= m12; s[F_]  += m12;
      real m13 = s[VF13_] * wfF_13; s[VF13_] -= m13; s[VF_] += m13;
      real m14 = s[NF13_] * wfF_13; s[NF13_] -= m14; s[NF_] += m14;
    }

    // aging + mortality (single Euler step per year) total_pop_before_step is summed across all 21 carriage families before this step, used for births in next step.
    real total_pop = 0;
    for (k in 0:20) for (i in 1:n) total_pop += s[k * n + i];

    int n_fam = 21;
    for (k in 0:(n_fam - 1)) {
      int base = k * n;
      vector[n] orig;
      for (i in 1:n) orig[i] = s[base + i];
      for (i in 1:n) {
        real out_rate = aging_rate[i] + mu_age[i];
        real new_i = orig[i] * (1.0 - out_rate);
        if (i > 1) new_i += orig[i - 1] * aging_rate[i - 1];
        s[base + i] = new_i;
      }
    }

    // vaccination at birth, split newborns 3 ways
    real births      = birth_rate * total_pop;
    real to_pcv13    = vacc_cov_pcv13 * births;
    real to_pcv7     = vacc_cov_pcv7  * births;
    real to_unvacc   = births - to_pcv7 - to_pcv13;
    s[1]           = s[1]           + to_unvacc;     // S<1y
    s[7  * n + 1]  = s[7  * n + 1]  + to_pcv7;       // S_7<1y
    s[14 * n + 1]  = s[14 * n + 1]  + to_pcv13;      // S_13<1y

    // reset incidence accumulators
    for (k in 21:23) for (i in 1:n) s[k * n + i] = 0.0;
    return s;
  }
}

// ============================================================================

data {
  int<lower=1> n_age;
  int<lower=1> n_year;
  int<lower=1> n_burn;
  int<lower=1, upper=n_year> n_fit_start_t;
  int<lower=1, upper=3> scenario;

  vector[n_age]              gamma_V;
  vector[n_age]              gamma_F;
  vector[n_age]              gamma_N;
  vector[n_age]              aging_rate;

  // explicit demographic inputs (from demographic_model.rds)
  vector[n_age]              mu_age;
  vector[n_year]             birth_rate_year;
  vector[n_age]              init_1999_fitted;

  matrix[n_age, n_age]       c_sym_1999;
  vector[n_age]              pop_1999;
  array[n_year] matrix[n_age, n_age] c_sym_annual;
  array[n_year] vector[n_age]        pop_annual;

  vector[24 * n_age]         seed_state;

  vector<lower=0, upper=1>[n_year] vacc_cov_pcv7_year;
  vector<lower=0, upper=1>[n_year] vacc_cov_pcv13_year;

  // joint-posterior MVN propagation on log scale
  // log_carr = [log_q_age(1..n_age), log_rr_V, log_rr_F, log_rr_N_pre]
  // log_ipd  = [log_delta_V_7, log_delta_F_7, log_omega_V_7, log_omega_F_7]
  // carr_L and ipd_L are lower-triangular Cholesky factors of each posterior covariance matrix, so multi_normal_cholesky(mean, L) samples log_x with
  // exact posterior mean, marginal SD and pairwise correlations. carr_lb/ub and ipd_lb/ub are per-element safety bounds at +-10 marginal
  // SDs -- statistically transparent under the MVN prior but they prevent exp() overflow in the ODE call during aggressive warmup leapfrogs.
  vector[n_age + 3]                   carr_mean;
  matrix[n_age + 3, n_age + 3]        carr_L;
  vector[n_age + 3]                   carr_lb;
  vector[n_age + 3]                   carr_ub;
  vector[4]                           ipd_mean;
  matrix[4, 4]                        ipd_L;
  vector[4]                           ipd_lb;
  vector[4]                           ipd_ub;

  // safety bounds for the (log-scale) estimated PCV13 parameters
  real                       log_delta_F_13_lb;
  real                       log_delta_F_13_ub;
  real                       log_omega_F_13_lb;
  real                       log_omega_F_13_ub;
  real                       log_rr_N_post_lb;
  real                       log_rr_N_post_ub;

  array[n_year, n_age, 3] int<lower=0> obs_ipd;
}

// =============================================================================

transformed data {
  int n_states  = 24 * n_age;
  real t0       = 0.0;
  array[1] real ts = {1.0};

  real ode_rel_tol      = 1e-3;
  real ode_abs_tol      = 1e-3;
  int  ode_max_num_steps = 10000;
}

// =============================================================================

parameters {
  // Joint-posterior propagated parameters with +-10 SD safety bounds
  // (bounds are per-element; MVN priors below give the actual shape).
  vector<lower=carr_lb, upper=carr_ub>[n_age + 3] log_carr;   // [log_q_age(1..n_age), log_rr_V, log_rr_F, log_rr_N_pre]
  vector<lower=ipd_lb,  upper=ipd_ub>[4]          log_ipd;    // [log_delta_V_7, log_delta_F_7, log_omega_V_7, log_omega_F_7]

  // Estimated PCV13 parameters (unchanged from previous version).
  // All three declared so the same binary serves all scenarios; only the
  // active one gets an informative prior in the model block.
  real<lower=log_delta_F_13_lb, upper=log_delta_F_13_ub> log_delta_F_13;
  real<lower=log_omega_F_13_lb, upper=log_omega_F_13_ub> log_omega_F_13;
  real<lower=log_rr_N_post_lb,  upper=log_rr_N_post_ub>  log_rr_N_post;
}

// =============================================================================

transformed parameters {
  // ---- Slice MVN vectors into the previously-named log-scale scalars ------
  vector[n_age] log_q_age    = head(log_carr, n_age);
  real          log_rr_V     = log_carr[n_age + 1];
  real          log_rr_F     = log_carr[n_age + 2];
  real          log_rr_N_pre = log_carr[n_age + 3];
  real          log_delta_V_7  = log_ipd[1];
  real          log_delta_F_7  = log_ipd[2];
  real          log_omega_V_7  = log_ipd[3];
  real          log_omega_F_7  = log_ipd[4];

  // PCV13 V efficacy / waning identical to PCV7
  real log_delta_V_13 = log_delta_V_7;
  real log_omega_V_13 = log_omega_V_7;

  // Per-scenario, ONE PCV13 F parameter is estimated and the OTHER TWO are
  // set equal to their PCV7 counterparts:
  // scenario 1 estimates delta_F_13 ; omega_F_13_used = omega_F_7 ; rr_N_post_used  = rr_N_pre
  // scenario 2 estimates omega_F_13 ; delta_F_13_used = delta_F_7 ; rr_N_post_used  = rr_N_pre
  // scenario 3 estimates rr_N_post  ; delta_F_13_used = delta_F_7 ; omega_F_13_used = omega_F_7
  real log_delta_F_13_used = (scenario == 1) ? log_delta_F_13 : log_delta_F_7;
  real log_omega_F_13_used = (scenario == 2) ? log_omega_F_13 : log_omega_F_7;
  real log_rr_N_post_used  = (scenario == 3) ? log_rr_N_post  : log_rr_N_pre;

  // natural-scale parameters
  vector[n_age] q_age    = exp(log_q_age);
  real rr_V              = exp(log_rr_V);
  real rr_F              = exp(log_rr_F);
  real rr_N_pre          = exp(log_rr_N_pre);

  // "raw estimated" natural-scale (only meaningful for the scenario where it is active; exposed so the informative-prior statement can attach to it)
  real rr_N_post         = exp(log_rr_N_post);
  real delta_F_13        = exp(log_delta_F_13);
  real omega_F_13        = exp(log_omega_F_13);

  // "used" natural-scale
  real delta_F_13_used   = exp(log_delta_F_13_used);
  real omega_F_13_used   = exp(log_omega_F_13_used);
  real rr_N_post_used    = exp(log_rr_N_post_used);

  // other natural-scale exports
  real delta_V_7         = exp(log_delta_V_7);
  real delta_F_7         = exp(log_delta_F_7);
  real omega_V_7         = exp(log_omega_V_7);
  real omega_F_7         = exp(log_omega_F_7);
  real delta_V_13        = exp(log_delta_V_13);
  real omega_V_13        = exp(log_omega_V_13);

  // output arrays
  array[n_year, n_age, 3] real pred_ipd;
  array[n_year, n_age, 3] real inc_year;
  array[n_year] vector[n_states] year_start_state;
  array[n_year] vector[n_states] year_end_state;
  matrix[n_age, 3] ccr_1999;
  matrix[n_age, 3] ccr_2010;

  vector[n_states] cur = seed_state;

  // in-Stan burn-in, dynamics with pre-PCV13 rr_N burn-in occurs entirely before 1999, so all PCV13 strata are empty and
  // the F_13 / N_post parameter values don't matter here; we still pass the *_used values for consistency with the fit loop.
  for (b in 1:n_burn) {
    array[1] vector[n_states] sol_b = ode_bdf_tol(pneumo_rhs, cur, t0, ts,
                                                  ode_rel_tol, ode_abs_tol,
                                                  ode_max_num_steps,
                                                  n_age, q_age,
                                                  gamma_V, gamma_F, gamma_N,
                                                  rr_V, rr_F, rr_N_pre,
                                                  delta_V_7, delta_F_7,
                                                  delta_V_13, delta_F_13_used,
                                                  c_sym_1999, pop_1999);
    cur = age_vaccinate(sol_b[1], n_age, aging_rate, mu_age,
                        0.0, 0.0,
                        omega_V_7, omega_F_7,
                        omega_V_13, omega_F_13_used,
                        birth_rate_year[1]);  // use 1999 birth rate during burn-in
  }

  // main fit loop: 1999..2019
  for (t in 1:n_year) {
    year_start_state[t] = cur;

    // period-specific rr_N: pre-PCV13 uses the propagated rr_N_pre
    // post-PCV13 uses rr_N_post_used (= rr_N_post if scenario 3, else rr_N_pre)
    real rr_N_t = (t <= 11) ? rr_N_pre : rr_N_post_used;
    array[1] vector[n_states] sol = ode_bdf_tol(pneumo_rhs, cur, t0, ts,
                                                ode_rel_tol, ode_abs_tol,
                                                ode_max_num_steps,
                                                n_age, q_age,
                                                gamma_V, gamma_F, gamma_N,
                                                rr_V, rr_F, rr_N_t,
                                                delta_V_7, delta_F_7,
                                                delta_V_13, delta_F_13_used,
                                                c_sym_annual[t], pop_annual[t]);
    vector[n_states] endv = sol[1];

    for (a in 1:n_age) {
      inc_year[t, a, 1] = endv[21 * n_age + a];
      inc_year[t, a, 2] = endv[22 * n_age + a];
      inc_year[t, a, 3] = endv[23 * n_age + a];
    }

    if (t == 1) {
      for (a in 1:n_age) {
        for (k in 1:3) {
          real inc_clip = fmax(inc_year[t, a, k], 1e-6);
          ccr_1999[a, k] = obs_ipd[t, a, k] / inc_clip;
        }
      }
    }

    // CCR_2010 = obs_2010 / inc_2009 (serotype-reshuffle for PCV13 era)
    if (t == 11) {
      for (a in 1:n_age) {
        for (k in 1:3) {
          real inc_clip = fmax(inc_year[t, a, k], 1e-6);
          ccr_2010[a, k] = obs_ipd[12, a, k] / inc_clip;
        }
      }
    }

    if (t <= 11) {
      for (a in 1:n_age) for (k in 1:3)
        pred_ipd[t, a, k] = inc_year[t, a, k] * ccr_1999[a, k];
    } else {
      for (a in 1:n_age) for (k in 1:3)
        pred_ipd[t, a, k] = inc_year[t, a, k] * ccr_2010[a, k];
    }

    year_end_state[t] = endv;
    cur = age_vaccinate(endv, n_age, aging_rate, mu_age,
                        vacc_cov_pcv7_year[t], vacc_cov_pcv13_year[t],
                        omega_V_7, omega_F_7,
                        omega_V_13, omega_F_13_used,
                        birth_rate_year[t]);
  }
}

// =============================================================================

model {
  // joint MVN propagation priors, these preserve every pairwise correlation from the source posteriors
  log_carr ~ multi_normal_cholesky(carr_mean, carr_L);
  log_ipd  ~ multi_normal_cholesky(ipd_mean,  ipd_L);

  // scenario-specific informative prior on the active PCV13 parameter prior is on the natural scale; the Jacobian for exp() is added.
  // the other two log_*_13 / log_rr_N_post parameters are still declared (so the same Stan binary can serve all three scenarios) but do not
  // affect the dynamics in their non-active scenario, the *_used variables in transformed parameters wire to the PCV7 propagated values.
  if (scenario == 1) {
    delta_F_13 ~ beta(3, 9);
    target += log_delta_F_13;
  } else if (scenario == 2) {
    omega_F_13 ~ gamma(3, 3);
    target += log_omega_F_13;
  } else if (scenario == 3) {
    rr_N_post ~ beta(7, 2);
    target += log_rr_N_post;
  }

  // Poisson likelihood over the PCV13 fit window (2010..2019)
  for (t in n_fit_start_t:n_year) {
    for (a in 1:n_age) {
      for (k in 1:3) {
        target += poisson_lpmf(obs_ipd[t, a, k] | fmax(pred_ipd[t, a, k], 0.0) + 1e-6);
      }
    }
  }
}

// =============================================================================

generated quantities {
  array[n_year, n_age, 3] real foi_year;
  array[n_year, n_age, 3] int  pred_ipd_rep;
  array[n_year, n_age, 3] real log_lik;

  for (t in 1:n_year) {
    vector[n_age] carrV;
    vector[n_age] carrF;
    vector[n_age] carrN;
    vector[n_states] yss = year_start_state[t];
    for (i in 1:n_age) {
      carrV[i] = yss[1*n_age + i] + yss[4*n_age + i] + yss[5*n_age + i]
                + yss[8*n_age + i] + yss[11*n_age + i] + yss[12*n_age + i]
                + yss[15*n_age + i] + yss[18*n_age + i] + yss[19*n_age + i];
      carrF[i] = yss[2*n_age + i] + yss[4*n_age + i] + yss[6*n_age + i]
                + yss[9*n_age + i] + yss[11*n_age + i] + yss[13*n_age + i]
                + yss[16*n_age + i] + yss[18*n_age + i] + yss[20*n_age + i];
      carrN[i] = yss[3*n_age + i] + yss[5*n_age + i] + yss[6*n_age + i]
                + yss[10*n_age + i] + yss[12*n_age + i] + yss[13*n_age + i]
                + yss[17*n_age + i] + yss[19*n_age + i] + yss[20*n_age + i];
    }
    vector[n_age] lamV = q_age .* (c_sym_annual[t] * (carrV ./ pop_annual[t]));
    vector[n_age] lamF = q_age .* (c_sym_annual[t] * (carrF ./ pop_annual[t]));
    vector[n_age] lamN = q_age .* (c_sym_annual[t] * (carrN ./ pop_annual[t]));
    for (a in 1:n_age) {
      foi_year[t, a, 1] = lamV[a];
      foi_year[t, a, 2] = lamF[a];
      foi_year[t, a, 3] = lamN[a];
      for (k in 1:3) {
        real lambda_k = fmax(pred_ipd[t, a, k], 0.0) + 1e-6;
        pred_ipd_rep[t, a, k] = poisson_rng(lambda_k);
        log_lik[t, a, k]      = poisson_lpmf(obs_ipd[t, a, k] | lambda_k);
      }
    }
  }
}
