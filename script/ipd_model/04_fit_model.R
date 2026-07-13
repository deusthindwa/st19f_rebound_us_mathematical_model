# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

#load the require packages
if (!require(pacman)){
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
}

pacman::p_load(
  char = c("tidyverse","remotes", "ggplot2", "dplyr","rio","tidyr","plyr","lubridate","reshape2","curl","patchwork", "posterior",
           "deSolve","adaptivetau","data.table","scales","readr","MASS","rootSolve", "labelVector","PropCIs", "bayesplot", "readr",
           "binom","coda","Rcpp","gmm","RcppArmadillo","devtools","lattice", "RColorBrewer", "lhs","png", "cmdstanr", "tibble",
           "readxl","viridis","zoo","lattice","latex2exp","ape","cowplot","gridExtra","grid","ggpubr","here","deSolve", "loo")
)

#choose "quick" or "production" run but the latter takes several days to complete
MCMC_MODE <- "quick"
SEED  <- 20260524
FIT_END_YEAR <- 2007 #toggle here to vary calibration window
ALL_YEARS  <- 1999:2009
VACC_COV_YEAR <- c(0.43, 0.99, 0.99, 0.99, rep(0.95, length(ALL_YEARS) - 4))
stopifnot(length(VACC_COV_YEAR) == length(ALL_YEARS))

#n_burn: in-Stan burn-in years using CURRENT draw's q_age and rr_*.
#20 years gives dynamics ample time to settle to per-draw steady state even for draws far from mid-prior values
N_BURN     <- 20

#run the following scripts prior
#source(here::here("script", "ipd_model", "01_data_prep.R"))
#source(here::here("script", "ipd_model", "02_pre_pcv_sim.R"))
source(here::here("01_data_prep.R"))
source(here::here("02_pre_pcv_sim.R"))
n_age <- fixed_data$n_age

#load explicit demographic calibration
if (!file.exists(here::here("results", "ipd_model", "demographic_model.rds"))) {
  #source(here::here("script", "ipd_model", "06_demographic_calibration.R"))
  source(here::here("06_demographic_calibration.R"))
}
demog <- 
  #readRDS(here::here("results", "ipd_model", "demographic_model.rds"))
  readRDS(here::here("demographic_model.rds"))

#calibrated mortality, calibrated 1999 init pop, birth rate at each year of the fit window (per capita per year).
mu_age_calib <- as.numeric(demog$mu_age) 
init_pop_1999  <- as.numeric(demog$init_1999_fitted)
birth_rate_year <- as.numeric(demog$birth_rate[as.character(ALL_YEARS)])
stopifnot(length(birth_rate_year) == length(ALL_YEARS), length(mu_age_calib) == n_age, length(init_pop_1999) == n_age)

#carriage-fit posterior 99th marginal quantiles
CARR_PRIOR_LO_PCTILE <- 0.005
CARR_PRIOR_HI_PCTILE <- 0.995

#building independent log-uniform priors from carriage_fit.rds
fit_carr_obj <- 
  #readRDS(here::here("results", "carriage_model", "carriage_fit.rds"))
  readRDS(here::here("carriage_fit.rds"))

#robust extraction: support a CmdStanMCMC/Fit object or an already-extracted draws_df / draws_array.
draws_carr <- tryCatch(
  posterior::as_draws_df(fit_carr_obj$draws()),
  error = function(e) posterior::as_draws_df(fit_carr_obj)
)

log_rho_cols <- paste0("log_rho_age[", seq_len(n_age), "]")
missing_rho <- setdiff(log_rho_cols, names(draws_carr))
if (length(missing_rho)) stop("carriage_fit.rds is missing: ", paste(missing_rho, collapse = ", "))

missing_eps <- setdiff(c("eps_V", "eps_F", "eps_N"), names(draws_carr))
if (length(missing_eps)) stop("carriage_fit.rds is missing: ", paste(missing_eps, collapse = ", "))

#require enough draws for stable quantiles (~16,000 e.g., 4000 samples per chain from 4 chains)
phi_mat <- as.matrix(draws_carr[, log_rho_cols])
phi_mat <- cbind(phi_mat, log_eps_V = log(draws_carr[["eps_V"]]), log_eps_F = log(draws_carr[["eps_F"]]), log_eps_N = log(draws_carr[["eps_N"]]))
colnames(phi_mat) <- c(paste0("log_q_age[", seq_len(n_age), "]"), "log_rr_V", "log_rr_F", "log_rr_N")
stopifnot(nrow(phi_mat) >= 200)

#marginal quantile bounds on log scale (one pair per parameter).
prior_q_lo_pct <- apply(phi_mat, 2, quantile, probs = CARR_PRIOR_LO_PCTILE, na.rm = TRUE)
prior_q_hi_pct <- apply(phi_mat, 2, quantile, probs = CARR_PRIOR_HI_PCTILE, na.rm = TRUE)

#marginal means (used only by R seed burn-in as point values).
phi_means <- colMeans(phi_mat)

log_q_lower <- as.numeric(prior_q_lo_pct[1:n_age])
log_q_upper <- as.numeric(prior_q_hi_pct[1:n_age])
log_rr_V_lb <- as.numeric(prior_q_lo_pct[n_age + 1])
log_rr_V_ub <- as.numeric(prior_q_hi_pct[n_age + 1])
log_rr_F_lb <- as.numeric(prior_q_lo_pct[n_age + 2])
log_rr_F_ub <- as.numeric(prior_q_hi_pct[n_age + 2])
log_rr_N_lb <- as.numeric(prior_q_lo_pct[n_age + 3])
log_rr_N_ub <- as.numeric(prior_q_hi_pct[n_age + 3])
stopifnot(all(log_q_upper > log_q_lower), log_rr_V_ub > log_rr_V_lb, log_rr_F_ub > log_rr_F_lb, log_rr_N_ub > log_rr_N_lb)

#use marginal posterior means as point values for lightweight R seed burn-in.
#only needs reasonable nuisance values to land near the pre-PCV steady state. 
#Stan re-equilibrates per draw inside its own burn-in.
q_mean_post <- exp(phi_means[1:n_age])
rr_V_post <- exp(phi_means[n_age + 1])
rr_F_post <- exp(phi_means[n_age + 2])
rr_N_post <- exp(phi_means[n_age + 3])
fixed_data$q_mean <- setNames(q_mean_post, fixed_data$age_labels)
fixed_data$rr_V <- rr_V_post
fixed_data$rr_F <- rr_F_post
fixed_data$rr_N <- rr_N_post

#lightweight R burn-in to produce a seed state
#use calibrated 1999 initial pop so seed state's carriage totals match demographic model's starting point.
seed_R <- simulate_pre_pcv_seed(n_burn_years = 40, verbose = FALSE, pop_99 = init_pop_1999)

#Stan state has 20 families per age (14 carriage + 6 incidence accumulators)
seed_state <- numeric(20 * n_age)
for (k in 0:6) {  #unvaccinated carriage families
  seed_state[(k * n_age + 1):((k + 1) * n_age)] <-
    seed_R[(k * n_age + 1):((k + 1) * n_age)]
}

#vaccinated carriage families (7..13) and 6 incidence accumulators (14..19) = 0
stopifnot(length(seed_state) == 20 * n_age)

#contact matrix (raw daily), we only pass raw daily contact matrix to Stan.
n_year <- length(ALL_YEARS)
N_FIT <- match(FIT_END_YEAR, ALL_YEARS)
contact_raw <- fixed_data$contact_raw

#observed IPD for years 1999..2009
obs <- fixed_data$obs_ipd_array
ipd_int <- array(0L, dim = c(n_year, n_age, 3))
for (t in seq_along(ALL_YEARS)) {
  yr <- as.character(ALL_YEARS[t])
  ipd_int[t, , ] <- as.integer(obs[yr, , ])
}

#bounds for delta and omega, passed to Stan as data
#parameter constraints in Stan set to 95 % prior credible interval; beta(13,7), gamma(9,75)
DELTA_BETA_A <- 13  
DELTA_BETA_B <- 7
OMEGA_GAMMA_SHAPE <- 4   
OMEGA_GAMMA_RATE <- 11
delta_prior_lo <- qbeta(CARR_PRIOR_LO_PCTILE,  DELTA_BETA_A, DELTA_BETA_B)
delta_prior_hi <- qbeta(CARR_PRIOR_HI_PCTILE,  DELTA_BETA_A, DELTA_BETA_B)
omega_prior_lo <- qgamma(CARR_PRIOR_LO_PCTILE, OMEGA_GAMMA_SHAPE, OMEGA_GAMMA_RATE)
omega_prior_hi <- qgamma(CARR_PRIOR_HI_PCTILE, OMEGA_GAMMA_SHAPE, OMEGA_GAMMA_RATE)

#Stan data list
stan_data <- list(
  n_age            = n_age,
  n_year           = n_year,
  n_fit            = N_FIT,
  n_burn           = N_BURN,
  gamma_V          = as.numeric(fixed_data$gamma_V),
  gamma_F          = as.numeric(fixed_data$gamma_F),
  gamma_N          = as.numeric(fixed_data$gamma_N),
  aging_rate       = as.numeric(fixed_data$aging_rate),
  contact_raw      = contact_raw,
  init_pop_1999    = init_pop_1999,
  birth_rate_year  = birth_rate_year,
  mu_age           = mu_age_calib,
  seed_state       = seed_state,
  log_q_lower      = log_q_lower,
  log_q_upper      = log_q_upper,
  log_rr_V_lb      = log_rr_V_lb, log_rr_V_ub = log_rr_V_ub,
  log_rr_F_lb      = log_rr_F_lb, log_rr_F_ub = log_rr_F_ub,
  log_rr_N_lb      = log_rr_N_lb, log_rr_N_ub = log_rr_N_ub,
  delta_prior_lo   = delta_prior_lo,
  delta_prior_hi   = delta_prior_hi,
  omega_prior_lo   = omega_prior_lo,
  omega_prior_hi   = omega_prior_hi,
  obs_ipd          = ipd_int,
  vacc_cov_year    = VACC_COV_YEAR)

#no init function but Stan defaults to Uniform(-2, 2) on unconstrained scale, spreading out points across each constrained parameter box. 
#this gives each of the 4 chains a different starting point and avoids all chains starting at midpoint

#compile and sample
spn_model <- 
  #cmdstanr::cmdstan_model(here::here('script', "ipd_model", "03_pneumo_pcv.stan"))
  cmdstanr::cmdstan_model(here::here("03_pneumo_pcv.stan"))

#run the STAN model HMC/NUTS
mcmc_args <- switch(
  MCMC_MODE,
  "quick" = list(
    iter_warmup = 500, 
    iter_sampling = 1000,
    chains = 4, 
    parallel_chains = 4,
    adapt_delta = 0.9, 
    max_treedepth = 12),
  
  "production" = list(
    iter_warmup = 1000, 
    iter_sampling = 4000,
    chains = 4, 
    parallel_chains = 4,
    adapt_delta = 0.95, 
    max_treedepth = 13),
  
  stop("MCMC_MODE must be 'quick' or 'production'")
)

ipd_fit <- do.call(spn_model$sample, c(
  list(data = stan_data, seed = SEED, refresh = 50),
  mcmc_args
))

#save fit files for post-processing.
#ipd_fit$save_object(file = here::here('results', "ipd_model", "ipd_fit.rds"))
ipd_fit$save_object(file = here::here("ipd_fit.rds"))

saveRDS(list(
  stan_data = stan_data,
  all_years = ALL_YEARS,
  n_fit = N_FIT,
  fit_end_year = FIT_END_YEAR,
  log_q_lower = log_q_lower,
  log_q_upper = log_q_upper,
  log_rr_V_lb = log_rr_V_lb, log_rr_V_ub = log_rr_V_ub,
  log_rr_F_lb = log_rr_F_lb, log_rr_F_ub = log_rr_F_ub,
  log_rr_N_lb = log_rr_N_lb, log_rr_N_ub = log_rr_N_ub,
  phi_means = phi_means,
  pctiles = c(lo = CARR_PRIOR_LO_PCTILE, hi = CARR_PRIOR_HI_PCTILE)), 
  #file = here::here('results', "ipd_model", "ipd_fit_inputs.rds"))
  file = here::here("ipd_fit_inputs.rds"))
