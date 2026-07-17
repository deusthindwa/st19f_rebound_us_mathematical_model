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

#04a_fit_scenario1.R
#scenario 1 (during PCV13, 2010-2019):
#estimate delta_F_13 ~ Beta(3, 9) on the natural scale (sampled on log scale; Jacobian applied in Stan)
#the other two PCV13 F parameters are set equal to their PCV7 counterparts inside the Stan model (in transformed parameters)
#omega_F_13 := omega_F_7 (propagated uncertainty, sampled from ipd_fit.rds posterior of log_omega_F_7)
#rr_N_13 := rr_N_7 (propagated uncertainty, sampled from carriage_fit.rds posterior of log_rr_N)

#source(here::here("script", "ipd_model_sensitivity", "04_fit_pcv13_common.R"))
source(here::here("04_fit_pcv13_common.R"))

stan_data <- c(common_stan_data, list(
  scenario = 1L,
  #active parameter, bounds wide so the informative Beta(3, 9) dominates.
  log_delta_F_13_lb = log(0.01),  
  log_delta_F_13_ub = log(0.99),
  #inactive declared parameters, bounds kept tight; the prior on these is uniform on the log scale 
  #do not enter the dynamics in scenario 1 (variables wire to log_omega_F_7, log_rr_N_pre instead).
  log_omega_F_13_lb = common_stan_data$log_omega_F_7_lb,  
  log_omega_F_13_ub = common_stan_data$log_omega_F_7_ub,
  log_rr_N_post_lb  = common_stan_data$log_rr_N_pre_lb,  
  log_rr_N_post_ub  = common_stan_data$log_rr_N_pre_ub
))

mod <- get_pcv13_model()
ipd_fit <- do.call(mod$sample, c(list(data = stan_data, seed = SEED_BASE + 1L, refresh = 50), mcmc_args))

ipd_fit$save_object(file = "ipd_fit_scenario1.rds")
saveRDS(list(stan_data = stan_data, all_years = ALL_YEARS, scenario = 1L), "ipd_fit_inputs_scenario1.rds")

#scenario 1 posterior (estimated parameter)"
print(ipd_fit$summary(c("delta_F_13")))

#scenario 1: PCV13 values that flow into the dynamics (= PCV7)
#omega_F_13_used = omega_F_7 and rr_N_post_used = rr_N_pre
print(ipd_fit$summary(c("delta_F_13_used", "omega_F_13_used", "rr_N_post_used")))

#scenario 1: propagated PCV7 parameters (carriage_fit + ipd_fit)
print(ipd_fit$summary(c("q_age", "rr_V", "rr_F", "rr_N_pre", "delta_V_7", "delta_F_7", "omega_V_7", "omega_F_7")))
