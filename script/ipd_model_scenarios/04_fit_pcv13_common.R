# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(tibble)
})

#choose "quick" or "production" run but the latter takes several days to complete
MCMC_MODE   <- "quick" #"quick" or "production"
SEED_BASE   <- 20260601
ALL_YEARS   <- 1999:2019 #21 years
N_BURN      <- 20 #in-Stan burn-in
N_FIT_START <- 12 #year index of 2010

#per-boundary vaccination coverage.
#1999 -> 2000 : 43 %  (2000 birth cohort)
#2000 -> 2001 : 98 %  (2001 birth cohorts) etc
# PCV7 ramps 43%, 98%, 98%, 98%, 95% for 2000, 2001 / 2002 / 2003+ cohorts,
# then stops at the 2009 -> 2010 boundary.
# PCV13 takes over at the 2009 -> 2010 boundary at the SAME steady-state
# coverage that PCV7 had reached (0.95) -- the switch is a vaccine swap,
# not a re-introduction, so there is no second ramp.
N <- length(ALL_YEARS)
VACC_COV_PCV7_YEAR  <- c(0.43, 0.99, 0.99, 0.99, rep(0.95, 6), rep(0.00, 11))
VACC_COV_PCV13_YEAR <- c(rep(0.00, 10), rep(0.95, 11))
stopifnot(length(VACC_COV_PCV7_YEAR)  == N)
stopifnot(length(VACC_COV_PCV13_YEAR) == N)
stopifnot(all(VACC_COV_PCV7_YEAR + VACC_COV_PCV13_YEAR <= 1))

#load other relevant scripts
#source(here::here("script", "ipd_model_scenarios", "01_data_prep.R"))
source(here::here("01_data_prep.R"))

#load demographic calibration
if (!file.exists(
  here::here("demographic_model.rds"))
  #here::here("results", "ipd_model_scenarios", "demographic_model.rds")
  ) {
  stop("06_demographic_calibration.rds not found -- expected user-supplied previous fit")
}

#source(here::here("script", "ipd_model_scenarios", "06_demographic_calibration.R"))
source(here::here("06_demographic_calibration.R"))

#demog <- readRDS(here::here("results", "ipd_model_scenarios", "demographic_model.rds"))
demog <- readRDS(here::here("demographic_model.rds"))

mu_age  <- as.numeric(demog$mu_age)
birth_rate_year  <- as.numeric(demog$birth_rate[as.character(ALL_YEARS)])
init_1999_fitted <- as.numeric(demog$init_1999_fitted)
model_pop_matrix <- demog$model_pop_matrix
stopifnot(length(mu_age)          == fixed_data$n_age)
stopifnot(length(birth_rate_year) == length(ALL_YEARS))
stopifnot(nrow(model_pop_matrix)  == length(ALL_YEARS))

#load posterior vectors
if (!file.exists(
  #here::here("results", "carriage_model", "carriage_fit.rds"))
  here::here("carriage_fit.rds"))
  ) {
  stop("carriage_fit.rds not found -- expected user-supplied previous fit.")
}

if (!file.exists(
  #here::here("results", "ipd_model", "ipd_fit.rds"))
  here::here("ipd_fit.rds"))
  ) {
  stop("ipd_fit.rds not found -- expected user-supplied previous fit.")
}

#carriage_fit <- readRDS(here::here("results", "carriage_model", "carriage_fit.rds"))
carriage_fit <- readRDS("carriage_fit.rds")

#ipd_fit <- readRDS(here::here("results", "ipd_model", "ipd_fit.rds"))
ipd_fit <- readRDS("ipd_fit.rds")

#extract a parameter's draws as a (n_draws x dim) numeric matrix.
.draws_mat <- function(fit, par) {
  arr <- fit$draws(par, format = "draws_matrix")
  as.matrix(arr)
}


#joint-posterior MVN propagation
#replace the previous eight independent log-uniform (marginal-quantile) priors with two multivariate-normal priors on the log scale. This
# preserves every pairwise correlation from the source posteriors
#carriage_fit -> log_carr = [log_q_age(1..n_age), log_rr_V, log_rr_F, log_rr_N_pre]
#ipd_fit -> log_ipd  = [log_delta_V_7, log_delta_F_7, log_omega_V_7, log_omega_F_7]
#carriage_fit stores eps_V/F/N on the natural scale, so we log-transform to align with Stan parameter names (log_rr_V/F/N_pre).
#ipd_fit already stores log_delta_V/F and log_omega_V/F on the log scale.

#carriage_fit joint posterior on log scale
carr_log_draws <- cbind(
  as.matrix(carriage_fit$draws("log_rho_age", format = "draws_matrix")),   # n_draws x n_age
  log(as.matrix(carriage_fit$draws("eps_V",   format = "draws_matrix"))),
  log(as.matrix(carriage_fit$draws("eps_F",   format = "draws_matrix"))),
  log(as.matrix(carriage_fit$draws("eps_N",   format = "draws_matrix")))
)
stopifnot(ncol(carr_log_draws) == fixed_data$n_age + 3L)
carr_mean <- as.numeric(colMeans(carr_log_draws))
carr_cov  <- cov(carr_log_draws)
carr_L    <- t(chol(carr_cov))                     # lower-triangular Cholesky

#ipd_fit joint posterior on log scale
ipd_log_draws <- cbind(
  as.matrix(ipd_fit$draws("log_delta_V",  format = "draws_matrix")),
  as.matrix(ipd_fit$draws("log_delta_F",  format = "draws_matrix")),
  as.matrix(ipd_fit$draws("log_omega_V",  format = "draws_matrix")),
  as.matrix(ipd_fit$draws("log_omega_F",  format = "draws_matrix"))
)
stopifnot(ncol(ipd_log_draws) == 4L)
ipd_mean <- as.numeric(colMeans(ipd_log_draws))
ipd_cov  <- cov(ipd_log_draws)
ipd_L    <- t(chol(ipd_cov))

#safety bounds at +-10 SD (per element)
#HMC leapfrog moves during early warmup pushes unconstrained log_carr / log_ipd far past reasonable values, causing exp() to overflow in the ODE call
#bounds at +-10 marginal SDs are ~15 orders of magnitude past the MVN prior mass but prevent overflow.
carr_sd <- sqrt(diag(carr_cov))
carr_lb <- carr_mean - 10 * carr_sd
carr_ub <- carr_mean + 10 * carr_sd
ipd_sd  <- sqrt(diag(ipd_cov))
ipd_lb  <- ipd_mean - 10 * ipd_sd
ipd_ub  <- ipd_mean + 10 * ipd_sd

#chain initialisation
#start log_carr and log_ipd at their MVN means (with tiny jitter to differentiate chains for R-hat). The scenario-active PCV13 parameter starts
#at `active_log_init` (log-mean of the informative Beta / Gamma prior); the two inactive PCV13 parameters start at the pre-PCV13 counterpart mean.
make_init_fun <- function(scenario, active_log_init, jitter_sd = 0.005) {
  # helpers to jitter and clip vectors/scalars safely inside declared bounds
  jit_v <- function(x) x + rnorm(length(x), 0, jitter_sd)
  clip_safe <- function(x, lb, ub, margin = 1e-3) {
    pmin(pmax(x, lb + margin), ub - margin)
  }
  function(chain_id) {
    # inactive-PCV13 defaults (natural log-scale means from the source fits)
    lF13 <- if (scenario == 1L) active_log_init else ipd_mean[2]                    # log_delta_F_7
    lO13 <- if (scenario == 2L) active_log_init else ipd_mean[4]                    # log_omega_F_7
    lRP  <- if (scenario == 3L) active_log_init else carr_mean[fixed_data$n_age + 3L]  # log_rr_N_pre
    list(
      log_carr       = jit_v(carr_mean),
      log_ipd        = jit_v(ipd_mean),
      log_delta_F_13 = clip_safe(jit_v(lF13), log(0.01), log(0.99)),
      log_omega_F_13 = clip_safe(jit_v(lO13), log(0.01), log(5.00)),
      log_rr_N_post  = clip_safe(jit_v(lRP),  log(0.01), log(0.99))
    )
  }
}

#R seed burn-in (still uses 1999 population + mean q/rr)
#source(here::here("script", "ipd_model_scenarios", "02_pre_pcv_sim.R"))
source("02_pre_pcv_sim.R")

seed_R <- simulate_pre_pcv_seed(n_burn_years = 40, verbose = FALSE)
n_age <- fixed_data$n_age

#pad to length 24 * n_age = 96; only the 28 unvacc carriage entries are non-zero.
seed_state <- numeric(24 * n_age)
for (k in 0:6) {
  seed_state[(k * n_age + 1):((k + 1) * n_age)] <-
    seed_R[(k * n_age + 1):((k + 1) * n_age)]
}

#replace age-1 S with the calibrated initial population so the
#explicit-demographics trajectory starts at the right pop level
#rescale each age group's compartments to match init_1999_fitted
for (i in 1:n_age) {
  total_a <- sum(sapply(0:6, function(k) seed_state[k * n_age + i]))
  if (total_a > 0) {
    sf <- init_1999_fitted[i] / total_a
    for (k in 0:6) seed_state[k * n_age + i] <- seed_state[k * n_age + i] * sf
  }
}

stopifnot(length(seed_state) == 24 * n_age)

#annual contact matrices and populations (from model_pop_matrix)
#use model-derived populations (from demographic calibration) so the contact-matrix symmetrisation and FOI denominator are consistent with
#the explicit-demographics dynamics inside the Stan model.
n_year <- length(ALL_YEARS)
c_sym_arr <- array(0, dim = c(n_year, n_age, n_age))
pop_arr <- array(0, dim = c(n_year, n_age))
for (t in seq_along(ALL_YEARS)) {
  yr        <- ALL_YEARS[t]
  pop_now   <- model_pop_matrix[as.character(yr), ]
  c_sym_d   <- symmetrise_contacts(fixed_data$contact_raw, pop_now)
  c_sym_arr[t, , ] <- c_sym_d * 365.25
  pop_arr[t, ]     <- pop_now
}
pop_1999 <- model_pop_matrix["1999", ]
c_sym_1999 <- symmetrise_contacts(fixed_data$contact_raw, pop_1999) * 365.25

#observed IPD cases
obs <- fixed_data$obs_ipd_array
ipd_int <- array(0L, dim = c(n_year, n_age, 3))
for (t in seq_along(ALL_YEARS)) {
  ipd_int[t, , ] <- as.integer(obs[as.character(ALL_YEARS[t]), , ])
}

#Stan data list (scenarios layer scenario id + estimated-param bounds) -
common_stan_data <- list(
  n_age            = n_age,
  n_year           = n_year,
  n_burn           = N_BURN,
  n_fit_start_t    = N_FIT_START,
  gamma_V          = as.numeric(fixed_data$gamma_V),
  gamma_F          = as.numeric(fixed_data$gamma_F),
  gamma_N          = as.numeric(fixed_data$gamma_N),
  aging_rate       = as.numeric(fixed_data$aging_rate),
  mu_age           = mu_age,
  birth_rate_year  = birth_rate_year,
  init_1999_fitted = init_1999_fitted,
  c_sym_1999       = c_sym_1999,
  pop_1999         = as.numeric(pop_1999),
  c_sym_annual     = c_sym_arr,
  pop_annual       = pop_arr,
  seed_state       = seed_state,
  vacc_cov_pcv7_year  = VACC_COV_PCV7_YEAR,
  vacc_cov_pcv13_year = VACC_COV_PCV13_YEAR,

  # Joint-posterior MVN inputs + per-element safety bounds
  carr_mean = carr_mean,
  carr_L    = carr_L,
  carr_lb   = carr_lb,
  carr_ub   = carr_ub,
  ipd_mean  = ipd_mean,
  ipd_L     = ipd_L,
  ipd_lb    = ipd_lb,
  ipd_ub    = ipd_ub,

  obs_ipd = ipd_int
)

mcmc_args <- switch(
  MCMC_MODE,
  "quick" = list(
    iter_warmup = 500,
    iter_sampling = 200,
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

get_pcv13_model <- function() {
  #cmdstanr::cmdstan_model(here::here("script", "ipd_model_scenarios", "03_pneumo_pcv13.stan"))
  cmdstanr::cmdstan_model(here::here("03_pneumo_pcv13.stan")) #force_recompile = TRUE
}
