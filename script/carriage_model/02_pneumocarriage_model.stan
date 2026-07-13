// by Deus Thindwa and Dan Weinberger
// age-structured mathematical model for pneumococcal transmission
// 05/23/2026

// ====================================================================

// pneumococcal_model.stan
//
// deterministic age-structured Susceptible-Infected-Susceptible
// (SIS) model of pneumococcal carriage transmission for three
// serotype groups (V = vaccine, F = 19F, N = non-vaccine) across
// four age groups (<1y, 1-4y, 5-17y, 18+y).
//
// seven free parameters with informative priors:
//
//   (a) log_rho_age[a]     4 age-specific log-baselines
//                          ~ Normal( rho_age_logmean[a], rho_age_logsd )
//
//                          Susceptibility is shared across serotypes:
//                          rho_V[i] = rho_F[i] = rho_N[i] = rho_age[i]
//
//   (b) eps_V, eps_F, eps_N     3 serotype-specific competition
//                               (relative risk of co-colonisation)
//                               in (0, 1)
//                               ~ Beta( eps_alpha, eps_beta )
//
//   prior centres for the age-baseline are drawn from pneumococcal
//   studies (e.g. Ojal 2017, Choi 2011/2012; Bottomley 2017), which report per-contact susceptibility of the order:
//       <1y    ~ 0.05      (log ~ -3.00)
//       1-4y   ~ 0.03      (log ~ -3.51)
//       5-17y  ~ 0.015     (log ~ -4.20)
//       18+y   ~ 0.005     (log ~ -5.30)
//
//   Beta(2,2) gives a unimodal weakly-informative prior on each
//   competition parameter, centred at 0.5 with sd ~0.224.  Use
//   eps_alpha = eps_beta = 1 if you want uniform U(0,1) instead.
//
// time units in years
// seven compartments per age group: S, V, F, N, VF, NV, NF.
// dual-carriage clearance: each component clears independently.
// aging modelled in proportions.
// ============================================================

functions {
  vector sis_ode(real t, vector y,
                 vector rho_V, vector rho_F, vector rho_N,
                 real   eps_V, real   eps_F, real   eps_N,
                 vector r_V,   vector r_F,   vector r_N,
                 vector l_age,
                 matrix C,
                 vector p_birth) {

    int n_age = num_elements(r_V);
    vector[7 * n_age] dy;

    vector[n_age] S;
    vector[n_age] V;
    vector[n_age] Fc;
    vector[n_age] Nc;
    vector[n_age] VF;
    vector[n_age] NV;
    vector[n_age] NF;
    vector[n_age] lam_V;
    vector[n_age] lam_F;
    vector[n_age] lam_N;

    for (i in 1:n_age) {
      S[i]  = y[(i - 1) * 7 + 1];
      V[i]  = y[(i - 1) * 7 + 2];
      Fc[i] = y[(i - 1) * 7 + 3];
      Nc[i] = y[(i - 1) * 7 + 4];
      VF[i] = y[(i - 1) * 7 + 5];
      NV[i] = y[(i - 1) * 7 + 6];
      NF[i] = y[(i - 1) * 7 + 7];
    }

    // force of infection (per year)
    for (i in 1:n_age) {
      real sV = 0;
      real sF = 0;
      real sN = 0;
      for (j in 1:n_age) {
        sV += C[i, j] * (V[j]  + NV[j] + VF[j]);
        sF += C[i, j] * (Fc[j] + NF[j] + VF[j]);
        sN += C[i, j] * (Nc[j] + NV[j] + NF[j]);
      }
      lam_V[i] = rho_V[i] * sV;
      lam_F[i] = rho_F[i] * sF;
      lam_N[i] = rho_N[i] * sN;
    }

    for (i in 1:n_age) {
      real rVi = r_V[i];
      real rFi = r_F[i];
      real rNi = r_N[i];
      real li  = l_age[i];

      // transmission
      real dS  = -(lam_V[i] + lam_F[i] + lam_N[i]) * S[i]
                 + rVi * V[i] + rFi * Fc[i] + rNi * Nc[i];

      real dV  =  lam_V[i] * S[i]
                 - rVi * V[i]
                 - eps_V * (lam_F[i] + lam_N[i]) * V[i]
                 + rFi * VF[i] + rNi * NV[i];

      real dFc =  lam_F[i] * S[i]
                 - rFi * Fc[i]
                 - eps_F * (lam_V[i] + lam_N[i]) * Fc[i]
                 + rVi * VF[i] + rNi * NF[i];

      real dNc =  lam_N[i] * S[i]
                 - rNi * Nc[i]
                 - eps_N * (lam_V[i] + lam_F[i]) * Nc[i]
                 + rVi * NV[i] + rFi * NF[i];

      real dVF =  eps_V * lam_F[i] * V[i]
                 + eps_F * lam_V[i] * Fc[i]
                 - (rVi + rFi) * VF[i];

      real dNV =  eps_V * lam_N[i] * V[i]
                 + eps_N * lam_V[i] * Nc[i]
                 - (rVi + rNi) * NV[i];

      real dNF =  eps_F * lam_N[i] * Fc[i]
                 + eps_N * lam_F[i] * Nc[i]
                 - (rFi + rNi) * NF[i];

      // aging
      if (i == 1) {
        dS  += li * (p_birth[1] - S[i]);
        dV  += li * (p_birth[2] - V[i]);
        dFc += li * (p_birth[3] - Fc[i]);
        dNc += li * (p_birth[4] - Nc[i]);
        dVF += li * (p_birth[5] - VF[i]);
        dNV += li * (p_birth[6] - NV[i]);
        dNF += li * (p_birth[7] - NF[i]);
      } else {
        dS  += li * (S[i - 1]  - S[i]);
        dV  += li * (V[i - 1]  - V[i]);
        dFc += li * (Fc[i - 1] - Fc[i]);
        dNc += li * (Nc[i - 1] - Nc[i]);
        dVF += li * (VF[i - 1] - VF[i]);
        dNV += li * (NV[i - 1] - NV[i]);
        dNF += li * (NF[i - 1] - NF[i]);
      }

      dy[(i - 1) * 7 + 1] = dS;
      dy[(i - 1) * 7 + 2] = dV;
      dy[(i - 1) * 7 + 3] = dFc;
      dy[(i - 1) * 7 + 4] = dNc;
      dy[(i - 1) * 7 + 5] = dVF;
      dy[(i - 1) * 7 + 6] = dNV;
      dy[(i - 1) * 7 + 7] = dNF;
    }

    return dy;
  }
}

// ============================================================

data {
  int<lower=1> n_age;                            // number of age groups (4)
  matrix[n_age, n_age] C;                        // contacts per year (population-symmetrised)
  vector<lower=0>[n_age] r_V;                    // clearance per year
  vector<lower=0>[n_age] r_F;
  vector<lower=0>[n_age] r_N;
  vector<lower=0>[n_age] l_age;                  // aging per year (= 1/d_years)
  vector[7] p_birth;                             // birth composition (all S)
  array[n_age] int<lower=0> N_total;             // sample size per age group
  array[n_age, 4] int<lower=0> y_obs;            // counts: S, V, F, N
  real t0;
  array[1] real ts;                              // equilibrium output time (years)


  // informative priors:
  //   log_rho_age[a]  ~ Normal( rho_age_logmean[a], rho_age_logsd )
  //   eps_X           ~ Beta( eps_alpha, eps_beta )
  vector[n_age] rho_age_logmean;                 // age-specific prior centres (log)
  real<lower=0> rho_age_logsd;                   // age prior sd (log scale)
  real<lower=0> eps_alpha;                       // Beta(a, b) shape 1  for eps_X
  real<lower=0> eps_beta;                        // Beta(a, b) shape 2  for eps_X

  // ODE tolerances
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1>  max_steps;
}

// ============================================================

transformed data {
  vector[7 * n_age] y0;
  for (i in 1:n_age) {
    y0[(i - 1) * 7 + 1] = 0.85;   // S
    y0[(i - 1) * 7 + 2] = 0.05;   // V
    y0[(i - 1) * 7 + 3] = 0.05;   // F
    y0[(i - 1) * 7 + 4] = 0.05;   // N
    y0[(i - 1) * 7 + 5] = 0;      // VF
    y0[(i - 1) * 7 + 6] = 0;      // NV
    y0[(i - 1) * 7 + 7] = 0;      // NF
  }
}

// ============================================================

parameters {
  // (a) 4 age-specific log-baselines carriage susceptibility
  vector<upper=2>[n_age] log_rho_age;            // exp(2) ~ 7.4 hard ceiling

  // (b) 3 serotype-specific competition parameters
  real<lower=0, upper=1> eps_V;
  real<lower=0, upper=1> eps_F;
  real<lower=0, upper=1> eps_N;
}

// ============================================================

transformed parameters {
  // age-only baseline susceptibility, shared across all three serotypes:
  vector<lower=0>[n_age] rho_age = exp(log_rho_age);

  vector<lower=0>[n_age] rho_V = rho_age;
  vector<lower=0>[n_age] rho_F = rho_age;
  vector<lower=0>[n_age] rho_N = rho_age;

  array[1] vector[7 * n_age] y_sol;
  array[n_age] vector[4] p_obs;
  array[n_age] vector[7] comp_eq;

  y_sol = ode_bdf_tol(sis_ode,
                      y0, t0, ts,
                      rel_tol, abs_tol, max_steps,
                      rho_V, rho_F, rho_N,
                      eps_V, eps_F, eps_N,
                      r_V,   r_F,   r_N,
                      l_age, C,    p_birth);

  for (i in 1:n_age) {
    real Si  = y_sol[1, (i - 1) * 7 + 1];
    real Vi  = y_sol[1, (i - 1) * 7 + 2];
    real Fi  = y_sol[1, (i - 1) * 7 + 3];
    real Ni  = y_sol[1, (i - 1) * 7 + 4];
    real VFi = y_sol[1, (i - 1) * 7 + 5];
    real NVi = y_sol[1, (i - 1) * 7 + 6];
    real NFi = y_sol[1, (i - 1) * 7 + 7];

    comp_eq[i][1] = Si;
    comp_eq[i][2] = Vi;
    comp_eq[i][3] = Fi;
    comp_eq[i][4] = Ni;
    comp_eq[i][5] = VFi;
    comp_eq[i][6] = NVi;
    comp_eq[i][7] = NFi;

    // split co-carriage incidence equally between component
    // single-carriage states before fitting.
    vector[4] p;
    p[1] = Si;
    p[2] = Vi + 0.5 * VFi + 0.5 * NVi;
    p[3] = Fi + 0.5 * VFi + 0.5 * NFi;
    p[4] = Ni + 0.5 * NVi + 0.5 * NFi;

    for (k in 1:4) p[k] = fmax(p[k], 1e-12);
    p = p / sum(p);
    p_obs[i] = p;
  }
}

// ============================================================

model {
  // informative priors
  // (a) 4 age-specific log-baselines
  for (a in 1:n_age) {
    log_rho_age[a] ~ normal(rho_age_logmean[a], rho_age_logsd);
  }
  // (b) 3 serotype-specific competition parameters in (0, 1)
  eps_V ~ beta(eps_alpha, eps_beta);
  eps_F ~ beta(eps_alpha, eps_beta);
  eps_N ~ beta(eps_alpha, eps_beta);

  // Likelihood
  for (i in 1:n_age) {
    target += multinomial_lpmf(y_obs[i] | p_obs[i]);
  }
}

// ============================================================

generated quantities {
  array[n_age, 4] int y_rep;        // posterior predictive replications
  vector[n_age] log_lik;            // per-age-group log-likelihood (LOO)
  vector[n_age] lambda_V_eq;        // equilibrium FOI (per year)
  vector[n_age] lambda_F_eq;
  vector[n_age] lambda_N_eq;

  for (i in 1:n_age) {
    y_rep[i]   = multinomial_rng(p_obs[i], N_total[i]);
    log_lik[i] = multinomial_lpmf(y_obs[i] | p_obs[i]);

    real sV = 0; real sF = 0; real sN = 0;
    for (j in 1:n_age) {
      sV += C[i, j] * (comp_eq[j][2] + comp_eq[j][6] + comp_eq[j][5]);
      sF += C[i, j] * (comp_eq[j][3] + comp_eq[j][7] + comp_eq[j][5]);
      sN += C[i, j] * (comp_eq[j][4] + comp_eq[j][6] + comp_eq[j][7]);
    }
    lambda_V_eq[i] = rho_V[i] * sV;
    lambda_F_eq[i] = rho_F[i] * sF;
    lambda_N_eq[i] = rho_N[i] * sN;
  }
}
