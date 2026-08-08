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

# 04b_fit_scenario2.R
#scenario 2 (during PCV13, 2010-2019):
#estimate omega_F_13 ~ Gamma(shape = 3, rate = 3) on the natural scale (sampled on log scale; Jacobian applied in Stan)
#the other two PCV13 F parameters are set equal to their PCV7 counterparts inside the Stan model (in transformed parameters)
#delta_F_13 := delta_F_7 (propagated joint MVN, ipd_fit posterior)
#rr_N_13 := rr_N_7 (propagated joint MVN, carriage_fit posterior)

#source(here::here("script", "ipd_model_scenarios", "04_fit_pcv13_common.R"))
source(here::here("04_fit_pcv13_common.R"))

stan_data <- c(common_stan_data, list(
  scenario = 2L,
  #all three PCV13 parameters use the same wide brackets. In scenario 2
  #only log_omega_F_13 is active (informative Gamma(3, 3) prior). The other
  #two are shadow declarations whose values never enter the dynamics
  #(the *_used variables in transformed parameters wire them to
  #log_delta_F_7 / log_rr_N_pre). Wide brackets are cheap because the
  #shadow posteriors are just uniform draws that have no downstream effect.
  log_delta_F_13_lb = log(0.01),  log_delta_F_13_ub = log(0.99),
  log_omega_F_13_lb = log(0.01),  log_omega_F_13_ub = log(5.00),
  log_rr_N_post_lb  = log(0.01),  log_rr_N_post_ub  = log(0.99)
))

#chain init: log_carr and log_ipd start at their MVN means; active
#omega_F_13 at the Gamma(3, 3) mean (= 1.0) on log scale.
init_fun <- make_init_fun(scenario = 2L, active_log_init = log(1.0))

mod <- get_pcv13_model()
ipd_fit_scen <- do.call(mod$sample, c(list(data = stan_data, seed = SEED_BASE + 2L, refresh = 50, init = init_fun), mcmc_args))

ipd_fit_scen$save_object(file = "ipd_fit_scenario2.rds")
saveRDS(list(stan_data = stan_data, all_years = ALL_YEARS, scenario = 2L), "ipd_fit_inputs_scenario2.rds")

#scenario 2 posterior (estimated parameter)
print(ipd_fit_scen$summary(c("omega_F_13")))

#scenario 2: PCV13 values that flow into the dynamics (= PCV7)
#delta_F_13_used = delta_F_7 and rr_N_post_used = rr_N_pre.
print(ipd_fit_scen$summary(c("delta_F_13_used", "omega_F_13_used", "rr_N_post_used")))

#scenario 2: propagated PCV7 parameters (carriage_fit + ipd_fit)
print(ipd_fit_scen$summary(c("q_age", "rr_V", "rr_F", "rr_N_pre", "delta_V_7", "delta_F_7", "omega_V_7", "omega_F_7")))
