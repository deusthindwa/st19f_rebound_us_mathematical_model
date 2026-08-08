
# Deterministic Age-Structured Susceptible-Infected-Susceptible (SIS) Pneumococcal Model

# PART A (Dynamics in Pre-PCV period)

- First, the pneumococcal model is fitted to English carriage data in Stan
- This pneumococcal model has 3 serotype groups (V = vaccine, F = 19F, N = non-vaccine) across 4 age groups (<1y, 1-4y, 5-17y, 18+y).
- Seven free parameters with informative priors:
  - (a) log_rho_age[a]: 4 age-specific log-baselines ~ Normal( rho_age_logmean[a], rho_age_logsd )
  - susceptibility is shared across serotypes: rho_V[i] = rho_F[i] = rho_N[i] = rho_age[i]
  - (b) eps_V, eps_F, eps_N: 3 serotype-specific competition (relative risk of co-colonisation) in (0, 1) ~ Beta( eps_alpha, eps_beta )
- Prior centres for age-baseline are drawn from pneumococcal studies (e.g. Ojal 2017, Choi 2011/2012; Bottomley 2017)
- Beta(2,2) gives a unimodal weakly-informative prior on each competition parameter, centred at 0.5 with sd ~0.224.
- Time units in years
- Seven compartments per age group: S, V, F, N, VF, NV, NF.
- Dual-carriage clearance: each component clears independently.
- Aging modelled as 1/(width of age group).


# PART B (Dynamics in pre-PCV & PCV7 periods)

## Introduction

- Similar dynamics as pre-PCV7 with addition of vaccination component
- R + Stan implementation of SIS model for pneumococcal carriage and IPD. 
- The Stan model uses `ode_bdf_tol` (backward differentiation) rather than `ode_rk45_tol` to manage stiffness of fit.
- 4 age groups (<1y, 1–4y, 5–17y, 18+y)
- 3 serotype classes (V = PCV7 vaccine type, F = serotype 19F, N = non-vaccine type), single and dual carriage.
- PCV7 introduced in **2000** in the US. Vaccination modelled **at birth** with coverage ramp.
- Birth cohort vaccination coverage: 2000 at 43%, 2001+ at 95%
- At every year boundary, fraction of newborns enters `S_7<1y` (vaccinated) and rest enter `S<1y` (unvaccinated). 
- The vaccine FOI reduction acts on infants for their full first year of life. 
- Older age groups acquire vaccinated individuals only by aging up, no additional vaccination is applied at any aging boundary. 
- 1999 is fully pre-vaccination; both 1999 and 2000 are part of the Poisson likelihood.

## The workflow
1. **The model simulates carriage incidence each year** (1999..2009) inside Stan
2. Model tracks **unvaccinated-stratum** and **vaccinated-stratum** incidence separately (`inc_year_unvacc`, `inc_year_vacc`)
3. **CCR_1999 is computed in Stan**: `CCR_1999[a, k] = obs_IPD_1999[a, k]/model_carriage_incidence_1999[a, k]`.
4. Because 1999 is pre-vaccination, total incidence equals unvaccinated incidence, so the CCR is well-defined.
5. **Predicted IPD** for year *t* and age *a* is simply total carriage incidence times the 1999 CCR
   - *V*: `pred_ipd = (inc_V_unvacc + inc_V_vacc) · CCR_V[1999, a]`
   - *F*: `pred_ipd = (inc_F_unvacc + inc_F_vacc) · CCR_F[1999, a]`
   - *N*: `pred_ipd = (inc_N_unvacc + inc_N_vacc) · CCR_N[1999, a]`  (no vaccine effect)
   - Vaccine effect on V, F operates through reduced FOI in vaccinated stratum inside the ODE 
   - In 1999 there are no vaccinated cohorts so `pred_ipd_1999 ≈ obs_ipd_1999` by construction
6. **Poisson likelihood** is summed over 9 years (1999..2007), 2008-2009 is OOS predictions.

## Parameters

### Estimated parameters (PCV7 carriage-blocking and waning)
1. All 4 estimated parameters live on **log scale** in Stan parameters block
2. (`log_delta_V`, `log_delta_F`, `log_omega_V`, `log_omega_F`).
3. The natural-scale variables are recovered as `exp(log_*)` in `transformed parameters` block.
4. Parameter constraints are tied to **prior 99% CrI**, computed as `qbeta(0.005/0.995, 13, 7)` and `qgamma(0.005/0.995, 2, 8)`
5. Constraints are passed to Stan as data, so the sampler only explores inside the prior's plausible region.

| Symbol | Name | Sampled as | Natural-scale | Prior | Role |
|---|---|---|---|---|---|
| δ_V,7 | `delta_V` | `log_delta_V` | prior 95 % CrI ≈ (0.48, 0.83) | **Beta(13, 7)** (mean 0.765) | **Estimated** |
| δ_F,7 | `delta_F` | `log_delta_F` | prior 95 % CrI ≈ (0.48, 0.83) | **Beta(13, 7)** (mean 0.65) | **Estimated** |
| ω_V,7 | `omega_V` | `log_omega_V` | prior 95 % CrI ≈ (0.03, 0.697) | **Gamma(2, 8)** (mean 0.25, duration ~4 yr) | **Estimated** |
| ω_F,7 | `omega_F` | `log_omega_F` | prior 95 % CrI ≈ (0.03, 0.697) | **Gamma(2, 8)** (mean 0.25, duration ~4 yr) | **Estimated** |

### Uncertainty-propagation parameters (carriage susceptibility and co-colonisation rr)
1. `q_age[1..4]`, `rr_V`, `rr_F`, `rr_N` are sampled independently on log scale,
2. Bounds are marginal 99th quantiles of carriage-fit posterior.
3. Aging rates (1, 1/4, 1/13, 1/59) are used

## Code workflow
1. (`06_demographic_calibration.R`) **Demographic calibration of mortality or net migration** .
   - An **explicit demographic model**: the four-compartment linear ODE.
   - Integrated (12 monthly Euler substeps per year) forward from the **given 1999 initial population**
   - `(<1y=258 310, 1-4y=1 033 240, 5-17y=3 497 789, 18+y=14 271 440)` **crude birth rates** (1999–2019)
   - 4 age-specific net mortality rates `mu_age = (mu_<1y, mu_1-4y, mu_5-17y, mu_18+y)` calibrated against observed pop counts.
   - Calibration is **every year of observed pop available** in `fixed_data$pop_matrix` intersected with `SIM_YEARS` (1999–2019) 
2. (`01_data_prep.R`) **Data preparation** keeps.
   - Unsymmetrised contact matrix
   - Observed IPD cases (1999-2009)
   - Fixed clearance rates `γ = 365.25 / duration_days`. 
   - The `q_lower/q_upper` and `rr_*_lower/rr_*_upper` 
3. (`02_pre_pcv_sim.R`) **R seed burn-in**
   - Runs 4-age × 7-compartment SIS ODE for 20 years with **posterior MEAN** of `q_age` / `rr_*` (from `carriage_fit.rds`) 
   - Runs fixed 1999 populations, with rescales pop, discrete aging  at each year boundary
   - Outputs a state vector that is already close to pre-PCV steady state
   - Stan re-equilibrates this state for each MCMC draw via its own short burn-in
4. (`03_pneumo_pcv.stan`) **Stan model**
   - **MVN priors on `q_age` and `rr_V/F/N`, includes cross-parameter correlation imposed
   - **Prior-CrI parameter bounds on δ/ω** formed in `transformed data` and used in the `parameters` block declarations
   - The prior shape (Beta(13, 7) on δ, Gamma(2, 8) on ω) applied in the model block with the exp() Jacobian
   - **Explicit demographic year-boundary update** (`age_demog_step`). 12 monthly Euler substeps integrate linear demog ODE for 1yr
   - Per substep, computes per-age pop from compartment sums, apply mortality, aging to all 14 carriage states uniformly per age 
   - Inject `birth_rate_year[t] × total_pop × dt` newborns into `S<1y` & `S_7<1y` by `vacc_cov_year[t]`. Pop evolves naturally
   - **In-Stan burn-in.** Uses `age_rescale_burnin` (no birth/death/aging) to hold pop at `init_pop_1999` while carriage equilibrates. 
   - Receives the R seed state and integrates `n_burn` more years (currently 20) using **current draw's** `q_age` and `rr_*`
   - The 1999 symmetrised contact matrix is built once in `transformed data` from `contact_raw` and `init_pop_1999`
   - **Fit loop** over years 1999..2009. At the start of each year the model recomputes per-age pop from compartments
   - Model then symmetrises `contact_raw` against it, so FOI normalisation & contact intensities are consistent with carriage state
   - Yr1 (1999) integrates 1yr of dynamics, resulting incidence is basis `ccr_1999[a, k]=obs_ipd_1999[a, k]/inc_1999_model[a, k]`
   - Year-boundary update is `age_demog_step`
   - **Predicted IPD.** `pred_ipd[t, a, k] = (inc_unvacc + inc_vacc) · CCR_1999[a, k]` for every year, age, and serotype.
   - **Generated quantities** include
     - `foi_year`, `carr_year_start`, `prev_year_start`, `carr_year_end`, `prev_year_end`, `model_pop_year_start`,
     - `model_pop_year_end`, `inc_year`, `pred_ipd`, `pred_ipd_rep`, `log_lik`, `ccr_1999`, `duration_V`, `duration_F`

5. (`04_fit_model.R`) **MCMC model initiation**.
   - Loads all relevant R scripts.
   - Builds Stan datalist with bound vectors, 1999 contact/pop block (for burn-in), per-year arrays, per-year vacc ramp, & prior bundles

5. (`05_post_process.R`) **Post-processing** produces
   - Chain mixing for 4 parameters on natural-scale and log-scale δ_V, δ_F, ω_V, ω_F
   - Pairwise posterior of δ_V, δ_F, ω_V, ω_F
   - Density overlays for 4 estimated parameters Beta(13, 7) / Gamma(2, 8)
   - Posterior CCRs in 1999 with 95 % CrI
   - Posterior FOI by year/age/serotype for **2000–2007**
   - Observed vs posterior-predicted IPD for 1999–2007
   - Population-weighted FOI 2000–2007
   - OOS posterior-predictive vs observed in 2008 and 2009
   - Modelled carriage prevalence (V, F, N) over time, **1999–2009** by age group with 95 % posterior CrI

6. (`07_pneumoipd_batch.sh`) **Batch run** on Yale's High-Performance Computing


# PART C (Dynamics in pre-PCV, PCV7 & PCV13 periods)

## Parameters
All sampled parameters live on the **log scale** in Stan parameters block 

| Symbol | Name (natural) | Sampled as | Natural-scale | Role |
|---|---|---|---|---|
| δ_V,7 | `delta_V` | `log_delta_V` | (0.806, 0.875) | Uncertainty propagation, MVN |
| δ_F,7 | `delta_F` | `log_delta_F` | (0.709, 0.791) | Uncertainty propagation, MVN |
| ω_V,7 | `omega_V` | `log_omega_V` | (0.244, 0.330) | Uncertainty propagation, MVN |
| ω_F,7 | `omega_F` | `log_omega_F` | (0.092, 0.958) | Uncertainty propagation, MVN |
| q_age[1] | `q_age[1]` | `log_q_age[1]` | (0.0125, 0.0174) | Uncertainty propagation, MVN |
| q_age[2] | `q_age[2]` | `log_q_age[2]` | (0.0150, 0.0185) | Uncertainty propagation, MVN |
| q_age[3] | `q_age[3]` | `log_q_age[3]` | (0.00864, 0.0108) | Uncertainty propagation, MVN |
| q_age[4] | `q_age[4]` | `log_q_age[4]` | (0.00517, 0.00635) | Uncertainty propagation, MVN |
| rr_V | `rr_V` | `log_rr_V` | (0.0712, 0.198) | Uncertainty propagation, MVN |
| rr_F | `rr_F` | `log_rr_F` | (0.352, 0.952)  | Uncertainty propagation, MVN |
| rr_N | `rr_N` | `log_rr_N` | (0.178, 0.524)  | Uncertainty propagation, MVN |

- The 4 PCV7 parameters are *uncertainty propagation* 
- The 4 `q_age` parameters are are *uncertainty propagation*
- The 3 `rr_*` parameters are *uncertainty-propagation*
- Sampled on the log scale from joint-posterior MVN for priors per MCMC iteration
- Posterior of eatimated PCV13 parms will reflect uncertainty in age susceptibility, in co-colonisation risks & vaccine effects. 
- `rr_V7`, `rr_F7`, `rr_N7` in the PCV7-vaccinated population share the same parameter as their unvaccinated counterparts

## Three PCV13 scenarios (1999-2019, but fitting to 2010-2019 only)
- (01_data_prep.R) fixed parameters, contact matrix, US populations 1999, observed IPD 1999-2019
- (02_pre_pcv_sim.R) R burn-in to produce a initial states of carriage pre-PCV
- (06_demographic_calibration.R) fits per-age net mortality & full population trajectory
- (03_pneumo_pcv13.stan) Stan model for 1999-2019 with the following 
  - Three vaccination strata (unvacc + PCV7 + PCV13) 
  - Explicit demographics (births, per-age mortality at year boundaries). 
  - Propagated parameters are sampled from posterior vectors of carriage_fit.rds and ipd_fit.rds via linear interpolation. 
  - Period-specific rr_N: rr_N_pre (propagated) drives 1999-2009; rr_N_post drives 2010-2019.
- (04_fit_pcv13_common.R) loads demographic_model.rds, carriage_fit.rds, ipd_fit.rds, and builds a shared Stan data list
  - (04a_fit_scenario1.R) scenario 1: estimate delta_F_13 ~ Beta(3, 9) -> fit_scenario1.rds
  - (04b_fit_scenario2.R) scenario 2: estimate omega_F_13 ~ Gamma(3, 3) -> fit_scenario2.rds
  - (04c_fit_scenario3.R) scenario 3: estimate rr_N_post  ~ Beta(7, 2) -> fit_scenario3.rds
- (05_post_process_pcv13.R) per-scenario diag, post-PCV13 carriage prev, counterfactual forward sim, cross-scenario WAIC

## Workflow detail
- The PCV13 model adds a **third vaccination stratum**: `S_13, V_13, F_13, N_13, VF_13, NV_13, NF_13`. 
- The Stan state grows from 80 to **96** states 
  - 21 carriage families + 3 total incidence accumulators per age × 4 ages). 
  - PCV7 vaccination stops in 2010 (`vacc_cov_pcv7_year = c(0.43, 0.95, …, 0.95, 0, 0, …)`)
  - PCV13 starts in 2010 (`vacc_cov_pcv13_year = c(0, …, 0, 0.95, 0.95, …)`). 
  - Existing PCV7 cohorts continue to wane and age into PCV13-era years.
- The first **CCRs.** `CCR_1999` computed inside Stan from year-1 incidence and applied to 1999–2009 predictions. 
- A second **CCRs.** `CCR_2010=obs_IPD_2010/inc_2009_model` computed in Stan at yr 11 & applied to 2010–2019 predictions.
- This reflects the serotype reshuffle from PCV7 to PCV13 (PCV13 covers 6 more serotypes).
- **Likelihood.** Summed over **2010–2019** (`n_fit_start_t = 12`). 1999–2009 dynamics do not enter the likelihood.
- **Three scenarios in one Stan file, three fits.** Each scenario estimates **exactly one** PCV13 F-related parameter; 
- The other two are **set equal to their PCV7 counterparts** (which are propagated from `ipd_fit`, `carriage_fit`):

| Scenario | Estimated (PCV13) | Prior | The other two are set to |
|---|---|---|---|
| 1 | δ_F^13 (PCV13 efficacy vs F) | Beta(3, 9) | ω_F^13 = ω_F^7,  rr_N^13 = rr_N^7 |
| 2 | ω_F^13 (PCV13 waning vs F)  | Gamma(3, 3) | δ_F^13 = δ_F^7,  rr_N^13 = rr_N^7 |
| 3 | rr_N^13 (post-PCV13 N rel. risk) | Beta(7, 2) | δ_F^13 = δ_F^7,  ω_F^13 = ω_F^7 |

All other parameters are **propagated by sampling from the previous posteriors**:
| Parameter | Source | Period |
|---|---|---|
| `log_q_age[1..4]`, `log_rr_V`, `log_rr_F` | `carriage_fit.rds` | pre + post (single sample per iter) |
| `log_rr_N` | `carriage_fit.rds` | **pre only** -> `rr_N_pre` |
| `log_delta_V_7`, `log_delta_V_13` (= V_7) | `ipd_fit.rds` | pre + post |
| `log_delta_F_7` | `ipd_fit.rds` | pre only |
| `log_omega_V_7`, `log_omega_V_13` (= V_7) | `ipd_fit.rds` | pre + post |
| `log_omega_F_7` | `ipd_fit.rds` | pre only |

-The three R scripts (`04a_/04b_/04c_`) are independent
- They all source `04_fit_pcv13_common.R` for shared Stan-data assembly, then layer on per-scenario bounds and `scenario` integer
- Independent `fit_scenarioN.rds` files saved, allowing 3 fits to be submitted as parallel HPC jobs (`07_pneumoipd_scenario_batch`)
