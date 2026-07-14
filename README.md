
# Deterministic Age-Structured Susceptible-Infected-Susceptible (SIS) Pneumococcal Model

# PART A (Dynamics during Pre-PCV7 period)

- In the first part of the work, the model pneumococcal model is fitted to the English carriage data in Stan
- SIS model of pneumococcal carriage transmission for three serotype groups (V = vaccine, F = 19F, N = non-vaccine) across four age groups (<1y, 1-4y, 5-17y, 18+y).
- Seven free parameters with informative priors:
  - (a) log_rho_age[a]: 4 age-specific log-baselines ~ Normal( rho_age_logmean[a], rho_age_logsd )
  - susceptibility is shared across serotypes: rho_V[i] = rho_F[i] = rho_N[i] = rho_age[i]
  - (b) eps_V, eps_F, eps_N: 3 serotype-specific competition (relative risk of co-colonisation) in (0, 1) ~ Beta( eps_alpha, eps_beta )
- Prior centres for the age-baseline are drawn from pneumococcal studies (e.g. Ojal 2017, Choi 2011/2012; Bottomley 2017)
- Beta(2,2) gives a unimodal weakly-informative prior on each competition parameter, centred at 0.5 with sd ~0.224.
- Time units in years
- Seven compartments per age group: S, V, F, N, VF, NV, NF.
- Dual-carriage clearance: each component clears independently.
- Aging modelled as 1/(width of age group).


# PART B (Dynamics during pre-PCV7 & PCV7 periods)

## Introduction

- Similar dynamics as pre-PCV7 with addition of vaccination component
- R + Stan implementation of the deterministic age-structured Susceptible-Infected-Susceptible model for pneumococcal carriage and IPD. 
- The Stan model uses `ode_bdf_tol` (backward differentiation) rather than `ode_rk45_tol` to manage stiffness of fit.
- Four age groups (<1y, 1–4y, 5–17y, 18+y)
- Three serotype classes (V = PCV7 vaccine type, F = serotype 19F, N = non-vaccine type), single and dual carriage.

- PCV7 introduced in **2000** in the US. Vaccination is modelled **at birth** with a four-step coverage ramp:
- Birth cohort, vaccination coverage: 2000 43%, 2001+ 95%

- At every year boundary, fraction of newborns enters `S_7<1y` (vaccinated) and rest enter `S<1y` (unvaccinated). 
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

| Symbol | Name | Sampled as | Natural-scale | Prior | Role |
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



# PART C (Dynamics during pre-PCV7, PCV7 & PCV13 periods)

## Parameters

All sampled parameters live on the **log scale** in Stan parameters block 

| Symbol | Name (natural) | Sampled as | Natural-scale | Role |
|---|---|---|---|---|
| δ_V,7 | `delta_V` | `log_delta_V` | (0.01, 0.99) safety; Beta(13, 7) prior | **Estimated**, mean 0.70 |
| δ_F,7 | `delta_F` | `log_delta_F` | (0.01, 0.99) safety; Beta(13, 7) prior | **Estimated**, mean 0.70 |
| ω_V,7 | `omega_V` | `log_omega_V` | (0.001, 5.0) safety; Gamma(9, 75) prior | **Estimated**, mean 0.12 (duration ~8.3 yr) |
| ω_F,7 | `omega_F` | `log_omega_F` | (0.001, 5.0) safety; Gamma(9, 75) prior | **Estimated**, mean 0.12 (duration ~8.3 yr) |
| q_age[1] | `q_age[1]` | `log_q_age[1]` | (0.0125, 0.0174) | Uncertainty propagation, log-uniform |
| q_age[2] | `q_age[2]` | `log_q_age[2]` | (0.0150, 0.0185) | Uncertainty propagation, log-uniform |
| q_age[3] | `q_age[3]` | `log_q_age[3]` | (0.00864, 0.0108) | Uncertainty propagation, log-uniform |
| q_age[4] | `q_age[4]` | `log_q_age[4]` | (0.00517, 0.00635) | Uncertainty propagation, log-uniform |
| rr_V | `rr_V` | `log_rr_V` | (0.0712, 0.198) | Uncertainty propagation, log-uniform |
| rr_F | `rr_F` | `log_rr_F` | (0.352, 0.952)  | Uncertainty propagation, log-uniform |
| rr_N | `rr_N` | `log_rr_N` | (0.178, 0.524)  | Uncertainty propagation, log-uniform |

The four PCV7 parameters are **estimated from the IPD data**. The four `q_age` entries and the three `rr_*` are *uncertainty-propagation* parameters — sampled on the log scale from log-uniform priors per MCMC iteration so the posterior of the four estimated parameters reflects parametric uncertainty in age-specific susceptibility and in the co-colonisation relative risks. `rr_V7`, `rr_F7`, `rr_N7` in the PCV7-vaccinated population share the same parameter as their unvaccinated counterparts (Table 1 makes them identical).

Aging rates (1, 1/4, 1/13, 1/59) are stored in `aging_rate` (named to avoid clashing with the δ symbol used for PCV7 carriage effectiveness).

## File order

### Workflow (1999-2019, three PCV13 scenarios; fitting on 2010-2019 only)

```
00_README.md              this file
01_data_prep.R            fixed parameters, contact matrix, US populations
                          1999-2019, observed IPD 1999-2019
02_pre_pcv_sim.R          lightweight R burn-in to produce a SEED state
06_demographic_calibration.R  fits per-age net mortality and the
                              full population trajectory; saves
                              demographic_model.rds (RUN FIRST)
03_pneumo_pcv13.stan      Stan model for 1999-2019 with three vaccination
                          strata (unvacc + PCV7 + PCV13) and EXPLICIT
                          DEMOGRAPHICS (births + per-age mortality at year
                          boundaries; no rescale). Propagated parameters
                          are sampled from posterior vectors of
                          carriage_fit.rds and ipd_fit.rds via linear
                          interpolation. Period-specific rr_N: rr_N_pre
                          (propagated) drives 1999-2009; rr_N_post drives
                          2010-2019.
04_fit_pcv13_common.R     Loads demographic_model.rds + carriage_fit.rds
                          + ipd_fit.rds, builds the shared Stan data list
04a_fit_scenario1.R       Scenario 1: estimate delta_F_13 ~ Beta(3, 9)
                                      -> fit_scenario1.rds
04b_fit_scenario2.R       Scenario 2: estimate omega_F_13 ~ Gamma(3, 3)
                                      -> fit_scenario2.rds
04c_fit_scenario3.R       Scenario 3: estimate rr_N_post  ~ Beta(7, 2)
                                      -> fit_scenario3.rds
05_post_process_pcv13.R   Per-scenario diagnostics + post-PCV13 carriage
                          prevalence plot + counterfactual forward sim
                          + cross-scenario WAIC comparison
```

### Required external inputs (user-supplied)

The driver expects **`carriage_fit.rds`** and **`ipd_fit.rds`** in the working directory. These are previous MCMC fits whose posterior draws contain:

- `carriage_fit.rds`: `log_q_age` (n_draws x 4), `log_rr_V`, `log_rr_F`, `log_rr_N`
- `ipd_fit.rds`: `log_delta_V_7`, `log_delta_F_7`, `log_omega_V_7`, `log_omega_F_7`

For each previous-fit parameter the driver computes the **marginal log-scale quantile bounds** `(quantile(log_x, LO), quantile(log_x, HI))` and passes them to Stan as `<lower, upper>` bounds on the corresponding Stan parameter. With no explicit prior in the model block, Stan's default uniform on the constrained log-scale parameter gives a **log-uniform prior on the natural scale** over the marginal CrI. Parameters are sampled **independently**; joint posterior correlations from the previous fits are not retained.

The quantile width is controlled by `CARR_PRIOR_LO/HI_PCTILE` (carriage params) and `IPD_PRIOR_LO/HI_PCTILE` (IPD params) in `04_fit_pcv13_common.R`:

| Setting | Meaning |
|---|---|
| `0.025 / 0.975` | 95% CrI (default) |
| `0.005 / 0.995` | 99% CrI (looser) |
| `0.050 / 0.950` | 90% CrI (stricter) |
| `0 / 1`         | full empirical (min, max) support |

### HPC submission pattern

Each scenario is independent and stateless once `01_data_prep.R` /
`02_pre_pcv_sim.R` have run, so the three scenarios can be submitted as
**three independent jobs**. A typical SLURM array job might look like:

```bash
#SBATCH --array=1-3
#SBATCH --cpus-per-task=4   # cmdstanr uses parallel_chains = 4
#SBATCH --mem=16G
#SBATCH --time=24:00:00

case $SLURM_ARRAY_TASK_ID in
  1) Rscript 04a_fit_scenario1.R ;;
  2) Rscript 04b_fit_scenario2.R ;;
  3) Rscript 04c_fit_scenario3.R ;;
esac
```

Each task writes its own `fit_scenarioN.rds` and `fit_inputs_scenarioN.rds`.
Run `Rscript 05_post_process_pcv13.R` once all three (or any subset) have
completed -- it skips scenarios whose fit files are missing.

## How to run

```r
setwd("pneumo_ipd")
# Step 0: required prerequisite -- demographic calibration
source("06_demographic_calibration.R")    # writes demographic_model.rds
# Step 1: copy your prior MCMC outputs into the folder
#   - carriage_fit.rds
#   - ipd_fit.rds
# Step 2: run the three PCV13 scenarios (independent; run in parallel on HPC)
source("04a_fit_scenario1.R")  # delta_F_13 ~ Beta(3, 9)
source("04b_fit_scenario2.R")  # omega_F_13 ~ Gamma(3, 3)
source("04c_fit_scenario3.R")  # rr_N_post  ~ Beta(7, 2)
source("05_post_process_pcv13.R")  # diagnostics + plots + WAIC comparison
```

`MCMC_MODE <- "quick"` at the top of `04_fit_pcv13_common.R` (500 warmup + 500 sampling) is enough to sanity-check the workflow; `MCMC_MODE <- "production"` (1000 + 1000, `adapt_delta = 0.95`) for a publication-quality posterior. `N_BURN` (default 20) controls the in-Stan burn-in years.

### Runtime / what to do if it's slow

A gradient evaluation for HMC integrates roughly `N_BURN + 11` ODE years for 11 parameters with autodiff, so the wall time depends mostly on:

| Lever | Default | Trade-off |
|---|---|---|
| `ode_rel_tol`, `ode_abs_tol` (in `03_pneumo_pcv13.stan`, `transformed data`) | `1e-4`, `1e-4` | Tighter tolerances ⇒ slower but more accurate ODE. The defaults give about 1 part in 10⁴ precision, well below the Poisson noise floor. |
| `N_BURN` (top of `04_fit_pcv13_common.R`) | 20 | More burn-in ⇒ slower per gradient. 20 years lets the dynamics fully settle for any draw inside the (now-narrow) parameter boxes. |
| `MCMC_MODE` (top of `04_fit_pcv13_common.R`) | `"quick"` (500 + 500) | Production halves step-size adaptation noise but doubles the iteration count. |

Indicative wall times for one scenario (4 chains in parallel, modern CPU):

* `MCMC_MODE = "quick"`, `N_BURN = 20`: ~1–2 hours.
* `MCMC_MODE = "production"`, `N_BURN = 20`: ~4–8 hours.

Running the three scenarios as a SLURM array gives total wall time ≈ longest single-scenario run. Slow warmup is normal — Stan's step-size and mass-matrix adaptation make the early iterations the most expensive.

## Dependencies

```r
install.packages(c("cmdstanr", "posterior", "bayesplot",
                   "deSolve", "ggplot2", "dplyr", "tidyr", "tibble"))
cmdstanr::install_cmdstan()
```

## Workflow detail

1. **Data preparation** (`01_data_prep.R`).
   Stores the US population table (1999–2019), the unsymmetrised contact matrix (Table 2 from the original brief), observed IPD case counts (1999–2019), the fixed clearance rates `γ = 365.25 / duration_days`, **bounds** (lower/upper) for `q_age[1..4]` and the three `rr_*`, and the per-scenario propagation ranges for the eight PCV-era parameters (δ_V^7, δ_V^13, δ_F^7, δ_F^13, ω_V^7, ω_V^13, ω_F^7, ω_F^13). The mean of `q_age` is stored too — used only by the lightweight R seed burn-in.

2. **R seed burn-in** (`02_pre_pcv_sim.R`).
   Runs the 4-age × 7-compartment SIS ODE for 40 years with **mean** `q_age` and `rr_*` and fixed 1999 populations, with discrete aging + rescaling at each year boundary. Outputs a near-steady-state seed vector that Stan re-equilibrates per MCMC draw via its in-Stan burn-in.

3. **Stan model** (`03_pneumo_pcv13.stan`).
   * **In-Stan burn-in.** Receives the R seed state and integrates `n_burn = 20` more years using the **current draw's** `q_age`, `rr_*`, and PCV-era parameters, with 1999 contact matrix and population and no vaccination. This places the dynamics at the per-draw steady state before year 1 (1999) begins.
   * **Fit loop** over years 1999..2019 (21 years). Year 1 incidence yields **`CCR_1999`**; year 11 incidence combined with year-12 observed IPD yields **`CCR_2010 = obs_ipd_2010 / inc_year_2009`** (per the brief, to handle the PCV7→PCV13 serotype reshuffle). Each year-boundary update applies waning of both PCV7 and PCV13 strata, aging, three-way splitting of newborns into unvacc / PCV7 / PCV13 by `vacc_cov_pcv7_year[t]` and `vacc_cov_pcv13_year[t]`, and rescaling to the next year's populations.
   * **Predicted IPD.**
     `pred_ipd[t, a, k] = inc_year[t, a, k] · CCR_1999[a, k]` for `t ≤ 11` (years 1999..2009),
     `pred_ipd[t, a, k] = inc_year[t, a, k] · CCR_2010[a, k]` for `t ≥ 12` (years 2010..2019).
     No post-ODE multiplier.
   * **Likelihood** sums Poisson terms only over `t = n_fit_start_t..n_year`, i.e. **2010..2019**. 1999–2009 dynamics are integrated for context, used to compute the two CCRs, and shown in obs-vs-pred plots, but do not contribute to the gradient.
   * **Generated quantities** export `foi_year`, `inc_year`, `pred_ipd`, `pred_ipd_rep`, `log_lik`, `ccr_1999`, `ccr_2010`, `duration_V_7`, `duration_V_13`, `duration_F_7`, `duration_F_13`.

4. **Shared assembly** (`04_fit_pcv13_common.R`). Builds the part of the Stan data list that is identical across scenarios: contact matrices, populations, target populations, seed state padded to length 96, two coverage vectors, MCMC settings (chains, warmup, sampling).

5. **Three scenario drivers** (`04a_fit_scenario1.R`, `04b_fit_scenario2.R`, `04c_fit_scenario3.R`).
   Each sources the common file, layers on scenario-specific bounds for `δ_F^13`, `ω_F^13`, and `rr_F`, sets a distinct `scenario` integer, runs `cmdstanr::sample()`, and writes its own `fit_scenarioN.rds` + `fit_inputs_scenarioN.rds`. The three drivers can be submitted as parallel HPC jobs because they share no mutable state.

6. **Post-processing** (`05_post_process_pcv13.R`). For each scenario whose fit file exists:
   * `scenarioN_plot_01_trace.png` — chain mixing for the estimated parameter (and its log-scale counterpart)
   * `scenarioN_plot_02_pairs.png` — pairwise posterior across estimated + key propagated parameters
   * `scenarioN_plot_03_prior_post.png` — log-uniform prior overlay vs posterior on the estimated parameter
   * `scenarioN_plot_04_obs_pred.png` — observed vs posterior-predicted IPD across 1999..2019
   * `scenarioN_ccr_2010_posterior.csv` — CCR_2010 posterior (median + 95 % CrI) by age × serotype
   * `scenarioN_parameter_summary.csv` — mean, median, 95 % CrI, R-hat, ESS

   If all three fits and the `loo` package are present, `scenario_waic_comparison.csv` and a printed `loo::loo_compare` table are written.

## Modelling notes / what changed from the previous iteration

* **No post-ODE progression-blocking factor.** `pred_ipd_X = (inc_X_unvacc + inc_X_vacc) · CCR_X_1999` for every serotype X. The vaccine effect on V/F is entirely the reduced FOI in the vaccinated stratum inside the ODE; there is no additional multiplier on IPD beyond carriage blocking. The ODE still tracks unvaccinated- and vaccinated-stratum incidence separately (6 accumulators per age, 24 of the 80 ODE states) — those are summed before applying the CCR.
* **Log-scale sampling with informative priors for δ and ω.** All four estimated parameters and the seven uncertainty-propagation parameters are sampled as `log_*` and recovered via `exp()` in `transformed parameters`. The priors are:
  - `delta_V`, `delta_F` ~ **Beta(14, 6)** on the natural scale (mean 0.70, 95 % CrI ≈ 0.50–0.87)
  - `omega_V`, `omega_F` ~ **Gamma(shape = 9, rate = 75)** on the natural scale (mean 0.12, 95 % CrI ≈ 0.06–0.21; mean duration of protection ≈ 8.3 years, 95 % CrI ≈ 4.8–18 years)
  - `q_age` and `rr_V`/`rr_F`/`rr_N` keep their log-uniform priors over the user-supplied bounds.
  A Jacobian `target += log_delta_V + log_delta_F + log_omega_V + log_omega_F` is added in the model block for the exp() transformations.
* **Vaccination at birth (per-year coverage).** Stan reads a `vacc_cov_year[n_year]` vector — one coverage value per year boundary. The driver currently passes `c(0.50, 0.90, 0.90, ..., 0.90)`, so the 2000 birth cohort is 50 % vaccinated and every birth cohort from 2001 onward is 90 %. The vaccine FOI reduction acts on infants throughout `<1y`, not only after they age into `1-4y`. Older age groups become partly vaccinated only by aging up. Changing this vector in R is the single point of edit for any per-year coverage scheme.
* **No init function.** Stan's default `Uniform(-2, 2)` initialisation on the unconstrained scale spreads the 4 chains across each parameter box, removing the "all chains start at midpoint" pinning that masked poor exploration in earlier runs.
* **`q_age` and `rr_V`/`rr_F`/`rr_N` are Stan parameters with log-uniform priors** over the per-parameter bounds — uncertainty propagation, not estimation. They are sampled jointly with `delta_*` and `omega_*` so the posterior of the four estimated parameters reflects parametric uncertainty in age-specific susceptibility and in co-colonisation relative risks.
* **In-Stan burn-in extended to `N_BURN = 20` years** so the dynamics fully settle to the per-draw steady state even when `q_age` and `rr_*` are far from their mid-prior values.
* **CCR_1999 is computed inside Stan** from each MCMC draw's 1999 carriage incidence (year 1999 has no vaccinated cohorts so total = unvaccinated incidence).
* **Carriage incidence drives IPD fitting:** in code terms, `inc_year_unvacc[t, a, k]` and `inc_year_vacc[t, a, k]` are accumulated by the ODE within each year, summed to give total incidence, and multiplied by `CCR_1999[a, k]` to give `pred_ipd[t, a, k]`, the Poisson rate compared to `obs_ipd[t, a, k]`.
* **Stiff ODE solver.** The Stan model uses `ode_bdf_tol` (backward differentiation) rather than `ode_rk45_tol`. Carriage dynamics with rapid clearance, dual-carriage coupling and a discrete year-boundary rescale form a mildly stiff system that RK45 handles poorly — BDF gives a substantial wall-time saving per gradient evaluation and avoids pre-warmup stalls.

## PCV13 extension (1999–2019, three scenarios)

The PCV13 model adds a **third vaccination stratum**: `S_13, V_13, F_13, N_13, VF_13, NV_13, NF_13`. The Stan state grows from 80 to **96** states (21 carriage families + 3 total incidence accumulators per age × 4 ages). PCV7 vaccination stops in 2010 (`vacc_cov_pcv7_year = c(0.50, 0.90, …, 0.90, 0, 0, …)`), and PCV13 starts in 2010 (`vacc_cov_pcv13_year = c(0, …, 0, 0.50, 0.90, …)`). Existing PCV7 cohorts continue to wane and age into PCV13-era years.

**CCRs.** `CCR_1999` is still computed inside Stan from year-1 incidence and applied to 1999–2009 predictions. A second case-carrier ratio `CCR_2010 = obs_IPD_2010 / inc_2009_model` is computed inside Stan at year 11 and applied to 2010–2019 predictions, reflecting the serotype reshuffle from PCV7 to PCV13 (PCV13 covers 6 more serotypes).

**Likelihood.** Summed only over **2010–2019** (`n_fit_start_t = 12`). 1999–2009 dynamics inform the obs-vs-pred plot and the CCR_2010 calculation but do not enter the likelihood.

**Three scenarios — one Stan file, three fits.** Each scenario estimates **exactly one** PCV13 F-related parameter; the other two are **set equal to their PCV7 counterparts** (which are propagated from `ipd_fit` / `carriage_fit`):

| Scenario | Estimated (PCV13) | Prior | The other two are set to |
|---|---|---|---|
| 1 | δ_F^13 (PCV13 efficacy vs F) | Beta(3, 9) | ω_F^13 = ω_F^7,  rr_N^13 = rr_N^7 |
| 2 | ω_F^13 (PCV13 waning vs F)  | Gamma(3, 3) | δ_F^13 = δ_F^7,  rr_N^13 = rr_N^7 |
| 3 | rr_N^13 (post-PCV13 N rel. risk) | Beta(7, 2) | δ_F^13 = δ_F^7,  ω_F^13 = ω_F^7 |

Inside `03_pneumo_pcv13.stan`'s `transformed parameters` block, the *_used* variables (`delta_F_13_used`, `omega_F_13_used`, `rr_N_post_used`) hold the value that actually flows into the dynamics. Each is conditionally assigned: if its scenario is active it equals the estimated Stan parameter; otherwise it equals the PCV7-propagated value. The Stan binary is the same for all three scenarios — only the `scenario` integer and the safety bounds differ between drivers.

All other parameters are **propagated by sampling from the previous posteriors**:

| Parameter | Source | Period |
|---|---|---|
| `log_q_age[1..4]`, `log_rr_V`, `log_rr_F` | `carriage_fit.rds` | pre + post (single sample per iter) |
| `log_rr_N` | `carriage_fit.rds` | **pre only** -> `rr_N_pre` |
| `log_delta_V_7`, `log_delta_V_13` (= V_7) | `ipd_fit.rds` | pre + post |
| `log_delta_F_7` | `ipd_fit.rds` | pre only |
| `log_omega_V_7`, `log_omega_V_13` (= V_7) | `ipd_fit.rds` | pre + post |
| `log_omega_F_7` | `ipd_fit.rds` | pre only |

In R: each previous-fit parameter's log-scale marginal `(LO, HI)` quantiles are computed and passed to Stan as `<lower, upper>` bounds. In Stan: each propagated parameter is declared `real<lower=lb, upper=ub>` with **no explicit prior** in the model block; Stan's default uniform on the constrained log-scale parameter gives a **log-uniform prior on the natural scale** over the marginal CrI. Parameters are sampled **independently** — joint correlations across propagated parameters are not retained, but their dependence with the scenario-estimated parameter is captured through the likelihood. Tightness of the propagation prior is controlled in R via `CARR_PRIOR_LO/HI_PCTILE` and `IPD_PRIOR_LO/HI_PCTILE` (95% / 99% / 90% CrI, or 0/1 for full empirical support).

The `scenario` integer in `stan_data` selects which estimated parameter gets the informative prior. The other two estimated parameters still exist (the Stan file compiles to one binary) but their log-scale safety bounds keep them tightly constrained so they cannot confound the inference of the scenario's target parameter.

The three drivers (`04a_/04b_/04c_`) are independent — they all source `04_fit_pcv13_common.R` for the shared Stan-data assembly, then layer on the per-scenario bounds and `scenario` integer, and save independent `fit_scenarioN.rds` files. This is what allows the three fits to be submitted as parallel HPC jobs.



