
# Pneumococcal IPD age-structured SIS model

## Introduction

- R + Stan implementation of the deterministic age-structured Susceptible-Infected-Susceptible model for pneumococcal carriage and IPD. 
- The Stan model uses `ode_bdf_tol` (backward differentiation) rather than `ode_rk45_tol` to manage stiffness of fit.
- Four age groups (<1y, 1–4y, 5–17y, 18+y), three serotype classes (V = PCV7 vaccine type, F = serotype 19F, N = non-vaccine type), single and dual carriage.

- PCV7 was introduced in **2000** in the US. Vaccination is modelled **at birth** with a four-step coverage ramp:
- Birth cohort, vaccination coverage: 2000 **43 %**, 2001 **98 %**, 2002, **98 %**, 2003 **98 %**, 2004+ **90 %**

- At every year boundary, that fraction of newborns enters `S_7<1y` (vaccinated) and the rest enter `S<1y` (unvaccinated). 
- The vaccine FOI reduction acts on infants for their full first year of life. 
- Older age groups acquire vaccinated individuals only by aging up, no additional vaccination is applied at any aging boundary. 
- 1999 is fully pre-vaccination; both 1999 and 2000 are part of the Poisson likelihood.

The flow is:
1. **The model simulates carriage incidence each year** (1999..2009) inside Stan, tracking **unvaccinated-stratum** and **vaccinated-stratum** incidence separately (`inc_year_unvacc`, `inc_year_vacc`)
2. **CCR_1999 is computed inside Stan**: `CCR_1999[a, k] = obs_IPD_1999[a, k] / model_carriage_incidence_1999[a, k]`. Because year 1999 is pre-vaccination, total incidence equals unvaccinated incidence, so the CCR is well-defined.
3. **Predicted IPD** for year *t* and age *a* is simply total carriage incidence times the 1999 CCR
   - *V*: `pred_ipd = (inc_V_unvacc + inc_V_vacc) · CCR_V[1999, a]`
   - *F*: `pred_ipd = (inc_F_unvacc + inc_F_vacc) · CCR_F[1999, a]`
   - *N*: `pred_ipd = (inc_N_unvacc + inc_N_vacc) · CCR_N[1999, a]`  (no vaccine effect)
   - Vaccine effect on V, F operates through reduced FOI in vaccinated stratum inside the ODE 
   - In 1999 there are no vaccinated cohorts so `pred_ipd_1999 ≈ obs_ipd_1999` by construction
4. **Poisson likelihood** is summed over all 11 years (1999..2009).

## Parameters

### Estimated parameters (PCV7 carriage-blocking and waning)

All four estimated parameters live on **log scale** in Stan parameters block (`log_delta_V`, `log_delta_F`, `log_omega_V`, `log_omega_F`). 
The natural-scale variables are recovered as `exp(log_*)` in `transformed parameters` block. 
The natural-scale parameter constraints are tied to **prior 99% CrI**, computed as `qbeta(0.005/0.995, 13, 7)` and `qgamma(0.005/0.995, 2, 8)` 
Constraints are passed to Stan as data, so the sampler can only explore inside the prior's plausible region.

| Symbol | Code name | Sampled as | Natural-scale bounds | Prior | Role |
|---|---|---|---|---|---|
| δ_V,7 | `delta_V` | `log_delta_V` | prior 95 % CrI ≈ (0.44, 0.83) | **Beta(13, 7)** (mean 0.765) | **Estimated** |
| δ_F,7 | `delta_F` | `log_delta_F` | prior 95 % CrI ≈ (0.48, 0.83) | **Beta(13, 7)** (mean 0.65) | **Estimated** |
| ω_V,7 | `omega_V` | `log_omega_V` | prior 95 % CrI ≈ (0.03, 0.697) | **Gamma(2, 8)** (mean 0.25, duration ~4 yr) | **Estimated** |
| ω_F,7 | `omega_F` | `log_omega_F` | prior 95 % CrI ≈ (0.03, 0.697) | **Gamma(2, 8)** (mean 0.25, duration ~4 yr) | **Estimated** |

### Uncertainty-propagation parameters (carriage susceptibility and co-colonisation rr)

`q_age[1..4]`, `rr_V`, `rr_F`, `rr_N` are sampled independently on log scale, with bounds are marginal 99th quantiles of carriage-fit posterior.

Aging rates (1, 1/4, 1/13, 1/59) are used


## Workflow detail

0. **Demographic calibration of mortality or net migration** (`06_demographic_calibration.R`).
   - An **explicit demographic model**: the four-compartment linear ODE.
   - Integrated (12 monthly Euler substeps per year) forward from the **given 1999 initial population** `(<1y=258 310, 1-4y=1 033 240, 5-17y=3 497 789, 18+y=14 271 440)` using **crude birth rates** (1999–2019). 
   - Four age-specific net mortality rates `mu_age = (mu_<1y, mu_1-4y, mu_5-17y, mu_18+y)` are calibrated by least squares (L-BFGS-B, relative-error SSE) against observed population counts.
   - The calibration window is **every year of observed pop available** in `fixed_data$pop_matrix` intersected with `SIM_YEARS` (1999–2019)
   - **Negative `mu_age` values.** The observed total population grows from 19.06 M (1999) to 28.55 M (2009), ~4.1 % per year, far above ~1.2 % per year supported by crude birth rate alone. 
   - The least-squares fit therefore lands on **negative `mu_age` values for all four age groups**, which is net demographic rate = biological mortality - net immigration. 
   - The magnitudes reflect a growing surveillance population.
   
   ```
   dN_1/dt = b(t)*Total - (alpha_1 + mu_1)*N_1
   dN_i/dt = alpha_{i-1}*N_{i-1} - (alpha_i + mu_i)*N_i   (i = 2, 3)
   dN_4/dt = alpha_3*N_3 - mu_4*N_4
   ```
 
1. **Data preparation** keeps (`01_data_prep.R`).
   - Population table (1999–2009, **calibration target only**)
   - Unsymmetrised contact matrix
   - Observed IPD cases (1999–2009)
   - Fixed clearance rates `γ = 365.25 / duration_days`. 
   - The `q_lower/q_upper` and `rr_*_lower/rr_*_upper` 

2. **R seed burn-in** (`02_pre_pcv_sim.R`).
   - Runs 4-age × 7-compartment SIS ODE for 20 years with **posterior MEAN** of `q_age` / `rr_*` (from `carriage_fit.rds`) 
   - Runs fixed 1999 populations, with discrete aging + rescaling at each year boundary. 
   - Outputs a state vector that is already close to pre-PCV steady state. 
   - Stan re-equilibrates this state for each MCMC draw via its own short burn-in.

3. **Stan model** (`03_pneumo_pcv.stan`).
   - **Independent log-uniform priors on `q_age` and `rr_V/F/N`, no cross-parameter correlation is imposed.
   - **Prior-CrI parameter bounds on δ/ω**. The log-scale equivalents are formed in `transformed data` and used in the `parameters` block declarations. 
   - The prior shape (Beta(13, 7) on δ, Gamma(2, 8) on ω) is still applied in the model block with the exp() Jacobian.
   - **Explicit demographic year-boundary update** (`age_demog_step`). 12 monthly Euler substeps integrate the linear demographic ODE for one year. 
   - Per substep: compute current per-age pop from compartment sums; apply mortality + aging to all 14 carriage families uniformly per age (carriage state preserved); 
   - Inject `birth_rate_year[t] × total_pop × dt` newborns into `S<1y` and `S_7<1y` split by `vacc_cov_year[t]`. pop evolves naturally.
   - **In-Stan burn-in.** Uses `age_rescale_burnin` (no birth/death/aging via rates) which holds the population at `init_pop_1999` while the carriage state equilibrates. 
   - Receives the R seed state and integrates `n_burn` more years (currently 20) using **current draw's** `q_age` and `rr_*`. The 1999 symmetrised contact matrix is built once in `transformed data` from `contact_raw` and `init_pop_1999`.
   - **Fit loop** over years 1999..2009. At the start of each year the model recomputes per-age pop from compartments and symmetrises `contact_raw` against it, so FOI normalisation and contact intensities are internally consistent with the carriage state. 
   - Year 1 (1999) integrates one year of dynamics, the resulting incidence becomes the basis for `ccr_1999[a, k] = obs_ipd_1999[a, k] / inc_1999_model[a, k]`. Year-boundary update is `age_demog_step`.
   - **Predicted IPD.** `pred_ipd[t, a, k] = (inc_unvacc + inc_vacc) · CCR_1999[a, k]` for every year, age, and serotype.
   - **Generated quantities** `foi_year`, `carr_year_start`, `prev_year_start`, `carr_year_end`, `prev_year_end`, `model_pop_year_start`, `model_pop_year_end`, `inc_year`, `pred_ipd`, `pred_ipd_rep`, `log_lik`, `ccr_1999`, `duration_V`, `duration_F`.

4. **MCMC model initiation** (`04_fit_model.R`).
   - Loads all relevant R scripts.
   - Builds the Stan data list with the bound vectors, the 1999 contact / population block (for burn-in), per-year arrays, the per-year vaccination ramp, and the prior bundles.

5. **Post-processing** produces (`05_post_process.R`).
   - Chain mixing for 4 parameters on natural-scale and log-scale δ_V, δ_F, ω_V, ω_F
   - Pairwise posterior of δ_V, δ_F, ω_V, ω_F
   - Density overlays for 4 estimated parameters Beta(13, 7) / Gamma(2, 8)
   - Posterior CCRs in 1999 with 95 % CrI
   - Posterior FOI by year/age/serotype for **2000–2007**
   - Observed vs posterior-predicted IPD for 1999–2009
   - Population-weighted FOI 2000–2007
   - Posterior-predictive vs observed in 2008 and 2009
   - Modelled carriage prevalence (V, F, N) over time, **1999–2009** by age group with 95 % posterior CrI

6. **Batch run** on Yale's High-Performance Computing (`07_pneumoipd_batch.sh`)
