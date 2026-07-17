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

#quantile width for the log-uniform "propagation prior" on pre-PCV13 parameter. 
#the Stan model gives each propagated parameter the prior log_x ~ Uniform( quantile(log_x_post, LO), quantile(log_x_post, HI) )
CARR_PRIOR_LO_PCTILE <- 0.005
CARR_PRIOR_HI_PCTILE <- 0.995
IPD_PRIOR_LO_PCTILE  <- 0.005
IPD_PRIOR_HI_PCTILE  <- 0.995

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
#source(here::here("script", "ipd_model_sensitivity", "01_data_prep.R"))
source(here::here("01_data_prep.R"))

#load demographic calibration
if (!file.exists(
  here::here("demographic_model.rds"))
  #here::here("results", "ipd_model_sensitivity", "demographic_model.rds")
  ) {
  stop("06_demographic_calibration.rds not found -- expected user-supplied previous fit")
}

#source(here::here("script", "ipd_model_sensitivity", "06_demographic_calibration.R"))
source(here::here("06_demographic_calibration.R"))

#demog <- readRDS(here::here("results", "ipd_model_sensitivity", "demographic_model.rds"))
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

#marginal log-scale quantile bounds for a vector of draws.
#Returns c(lb, ub); also enforces lb < ub (small jitter if degenerate).
.q_bounds <- function(x, lo, hi) {
  #stopifnot(length(x) >= 2L, lo >= 0, hi <= 1, lo < hi)
  b <- if (lo == 0 && hi == 1) c(min(x), max(x))
       else as.numeric(quantile(x, probs = c(lo, hi), na.rm = TRUE))
  if (diff(b) <= 0) b[2] <- b[1] + 1e-6
  b
}

#marginal log-scale quantile bounds from carriage_fit posterior
.log_q_age_post <- .draws_mat(carriage_fit, "log_rho_age")  #n_draws x n_age
stopifnot(ncol(.log_q_age_post) == fixed_data$n_age)
log_q_age_bnds <- apply(.log_q_age_post, 2, .q_bounds, lo = CARR_PRIOR_LO_PCTILE, hi = CARR_PRIOR_HI_PCTILE)  #2 x n_age
log_q_age_lb <- as.numeric(log_q_age_bnds[1, ])
log_q_age_ub <- as.numeric(log_q_age_bnds[2, ])

log_rr_V_bnds <- log(.q_bounds(as.vector(.draws_mat(carriage_fit, "eps_V")), CARR_PRIOR_LO_PCTILE, CARR_PRIOR_HI_PCTILE))
log_rr_F_bnds <- log(.q_bounds(as.vector(.draws_mat(carriage_fit, "eps_F")), CARR_PRIOR_LO_PCTILE, CARR_PRIOR_HI_PCTILE))
log_rr_N_pre_bnds <- log(.q_bounds(as.vector(.draws_mat(carriage_fit, "eps_N")), CARR_PRIOR_LO_PCTILE, CARR_PRIOR_HI_PCTILE))

#marginal log-scale quantile bounds from ipd_fit posterior
log_delta_V_7_bnds <- .q_bounds(as.vector(.draws_mat(ipd_fit, "log_delta_V")), IPD_PRIOR_LO_PCTILE, IPD_PRIOR_HI_PCTILE)
log_delta_F_7_bnds <- .q_bounds(as.vector(.draws_mat(ipd_fit, "log_delta_F")), IPD_PRIOR_LO_PCTILE, IPD_PRIOR_HI_PCTILE)
log_omega_V_7_bnds <- .q_bounds(as.vector(.draws_mat(ipd_fit, "log_omega_V")), IPD_PRIOR_LO_PCTILE, IPD_PRIOR_HI_PCTILE)
log_omega_F_7_bnds <- .q_bounds(as.vector(.draws_mat(ipd_fit, "log_omega_F")), IPD_PRIOR_LO_PCTILE, IPD_PRIOR_HI_PCTILE)

#R seed burn-in (still uses 1999 population + mean q/rr)
#source(here::here("script", "ipd_model_sensitivity", "02_pre_pcv_sim.R"))
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
  log_q_age_lb      = log_q_age_lb, #marginal log-scale quantile bounds (carriage_fit)
  log_q_age_ub      = log_q_age_ub,
  log_rr_V_lb       = log_rr_V_bnds[1],
  log_rr_V_ub       = log_rr_V_bnds[2],
  log_rr_F_lb       = log_rr_F_bnds[1],
  log_rr_F_ub       = log_rr_F_bnds[2],
  log_rr_N_pre_lb   = log_rr_N_pre_bnds[1],
  log_rr_N_pre_ub   = log_rr_N_pre_bnds[2],
  log_delta_V_7_lb  = log_delta_V_7_bnds[1], #marginal log-scale quantile bounds (ipd_fit)
  log_delta_V_7_ub  = log_delta_V_7_bnds[2],
  log_delta_F_7_lb  = log_delta_F_7_bnds[1],
  log_delta_F_7_ub  = log_delta_F_7_bnds[2],
  log_omega_V_7_lb  = log_omega_V_7_bnds[1],
  log_omega_V_7_ub  = log_omega_V_7_bnds[2],
  log_omega_F_7_lb  = log_omega_F_7_bnds[1],
  log_omega_F_7_ub  = log_omega_F_7_bnds[2],
  obs_ipd           = ipd_int
)

mcmc_args <- switch(
  MCMC_MODE,
  "quick" = list(
    iter_warmup = 500, 
    iter_sampling = 500,
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
  #cmdstanr::cmdstan_model(here::here("script", "ipd_model_sensitivity", "03_pneumo_pcv13.stan"))
  cmdstanr::cmdstan_model(here::here("03_pneumo_pcv13.stan"))
}
