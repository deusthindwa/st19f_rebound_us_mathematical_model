# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

# Extend the model's out-of-sample IPD projection past the Stan-model horizon (last year in ALL_YEARS, currently 2009) up to PROJ_END_YEAR (2019), without re-running Stan.
#
# Approach:
# 1. extract every posterior draw of the model parameters and of `year_end_state[n_year]` (state at the end of 2009).
# 2. re-implement Stan's dynamic model (within-year SIS ODE) and `age_vaccinate_demog` (year-boundary birth + death + aging + vaccine injection + waning) in R.
# 3. for each posterior draw of the 500 draws:
#   -apply the boundary step that Stan performs at the end of its last year (birth_rate_year[n_year], VACC_COV_YEAR[n_year]) to move from end-of-2009 to start-of-2010.
#   -iterate ODE + boundary-step year-by-year through PROJ_END_YEAR.
#   -multiply the resulting per-year, per-age, per-serotype incidence by the same draw's `ccr_1999` to get `pred_ipd[year, age, serotype]`.
# 4. aggregate across draws to get posterior-predictive intervals and combine with Stan's own in-sample + short-Out-of-Sample (2008-2009) `pred_ipd_rep` for a seamless 1999..PROJ_END_YEAR plot.
#
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(deSolve)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

source(here::here("script", "ipd_model", "01_data_prep.R"))

fit <- readRDS(here::here("results", "ipd_model", "ipd_fit.rds"))
inputs <- readRDS(here::here("results", "ipd_model", "ipd_fit_inputs.rds"))
all_years <- inputs$all_years
n_year_stan <- length(all_years)
demog <- readRDS(here::here("results", "ipd_model", "demographic_model.rds"))

#configurations
PROJ_END_YEAR  <- 2019
PROJ_YEARS     <- (max(all_years) + 1):PROJ_END_YEAR #2010..2019
VACC_COV_PROJ  <- rep(0.95, length(PROJ_YEARS)) #constant post-2009
N_PROJ_DRAWS   <- 500L
ODE_RTOL       <- 1e-4
ODE_ATOL       <- 1e-4
N_SUBSTEPS     <- 12L #matches Stan

#fixed data pulled from 01_data_prep.R and demographic_model.rds
n_age       <- fixed_data$n_age
n_fam       <- 20  #14 carriage + 6 inc
n_states    <- n_fam * n_age
age_labels  <- fixed_data$age_labels
serotypes   <- c("V", "F", "N")
contact_raw <- fixed_data$contact_raw
aging_rate  <- as.numeric(fixed_data$aging_rate)
gamma_V     <- as.numeric(fixed_data$gamma_V)
gamma_F     <- as.numeric(fixed_data$gamma_F)
gamma_N     <- as.numeric(fixed_data$gamma_N)
mu_age      <- as.numeric(demog$mu_age)
birth_rate  <- demog$birth_rate 

#model dynamics in R (mirror of Stan)
#within-year ODE
pneumo_rhs_R <- function(t, y, parms) {
  n   <- parms$n_age
  seg <- function(k) y[(k * n + 1):((k + 1) * n)]
  S <- seg(0);  V  <- seg(1);  Fc <- seg(2);  Nc <- seg(3)
  VF <- seg(4); NV <- seg(5);  NF <- seg(6)
  S7 <- seg(7); V7 <- seg(8);  F7 <- seg(9);  N7 <- seg(10)
  VF7 <- seg(11); NV7 <- seg(12); NF7 <- seg(13)
  
  carrV <- V  + VF  + NV  + V7  + VF7 + NV7
  carrF <- Fc + VF  + NF  + F7  + VF7 + NF7
  carrN <- Nc + NV  + NF  + N7  + NV7 + NF7
  
  lamV <- parms$q_age * as.vector(parms$c_sym %*% (carrV / parms$pop_j))
  lamF <- parms$q_age * as.vector(parms$c_sym %*% (carrF / parms$pop_j))
  lamN <- parms$q_age * as.vector(parms$c_sym %*% (carrN / parms$pop_j))
  
  lamV7 <- (1 - parms$delta_V) * lamV
  lamF7 <- (1 - parms$delta_F) * lamF
  lamN7 <- lamN
  
  gV <- parms$gamma_V; gF <- parms$gamma_F; gN <- parms$gamma_N
  rV <- parms$rr_V;   rF <- parms$rr_F;   rN <- parms$rr_N
  #rr's in the vaccinated stratum are the same (per Stan)
  rV7 <- rV; rF7 <- rF; rN7 <- rN
  
  dS  <- -S  * (lamV + lamF + lamN) + gV * V  + gF * Fc + gN * Nc
  dV  <-  S  * lamV - V  * (gV + rV * (lamF + lamN)) + gF * VF + gN * NV
  dF  <-  S  * lamF - Fc * (gF + rF * (lamV + lamN)) + gV * VF + gN * NF
  dN  <-  S  * lamN - Nc * (gN + rN * (lamV + lamF)) + gV * NV + gF * NF
  dVF <-  rV * V  * lamF + rF * Fc * lamV - VF * (gV + gF)
  dNV <-  rV * V  * lamN + rN * Nc * lamV - NV * (gV + gN)
  dNF <-  rN * Nc * lamF + rF * Fc * lamN - NF * (gN + gF)
  
  dS7  <- -S7 * (lamV7 + lamF7 + lamN7) + gV * V7  + gF * F7 + gN * N7
  dV7  <-  S7 * lamV7 - V7 * (gV + rV7 * (lamF7 + lamN7)) + gF * VF7 + gN * NV7
  dF7  <-  S7 * lamF7 - F7 * (gF + rF7 * (lamV7 + lamN7)) + gV * VF7 + gN * NF7
  dN7  <-  S7 * lamN7 - N7 * (gN + rN7 * (lamV7 + lamF7)) + gV * NV7 + gF * NF7
  dVF7 <-  rV7 * V7 * lamF7 + rF7 * F7 * lamV7 - VF7 * (gV + gF)
  dNV7 <-  rV7 * V7 * lamN7 + rN7 * N7 * lamV7 - NV7 * (gV + gN)
  dNF7 <-  rN7 * N7 * lamF7 + rF7 * F7 * lamN7 - NF7 * (gN + gF)
  
  #six incidence accumulators (unvaccinated- and vaccinated-stratum V/F/N)
  d_incV_u <-  S  * lamV  + rF  * Fc * lamV  + rN  * Nc * lamV
  d_incF_u <-  S  * lamF  + rV  * V  * lamF  + rN  * Nc * lamF
  d_incN_u <-  S  * lamN  + rV  * V  * lamN  + rF  * Fc * lamN
  d_incV_v <- S7 * lamV7 + rF7 * F7 * lamV7 + rN7 * N7 * lamV7
  d_incF_v <- S7 * lamF7 + rV7 * V7 * lamF7 + rN7 * N7 * lamF7
  d_incN_v <- S7 * lamN7 + rV7 * V7 * lamN7 + rF7 * F7 * lamN7
  
  dy <- numeric(20 * n)
  put <- function(k, v) dy[(k * n + 1):((k + 1) * n)] <<- v
  put(0,  dS);   put(1,  dV);   put(2,  dF);   put(3,  dN)
  put(4,  dVF);  put(5,  dNV);  put(6,  dNF)
  put(7,  dS7);  put(8,  dV7);  put(9,  dF7);  put(10, dN7)
  put(11, dVF7); put(12, dNV7); put(13, dNF7)
  put(14, d_incV_u); put(15, d_incF_u); put(16, d_incN_u)
  put(17, d_incV_v); put(18, d_incF_v); put(19, d_incN_v)
  list(dy)
}


#year-boundary step, matching `age_vaccinate_demog` in Stan (Euler over
#n_substeps monthly ticks: waning -> mortality+aging -> births).
age_vaccinate_demog_R <- function(y, n_age, aging_rate, mu_age,
                                  vacc_cov, omega_V, omega_F, birth_rate_y,
                                  n_substeps = N_SUBSTEPS) {
  n <- n_age
  s <- y
  wfV <- 1 - exp(-omega_V)
  wfF <- 1 - exp(-omega_F)
  dt  <- 1 / n_substeps
  
  #waning: vaccinated compartments leak to their unvaccinated counterparts
  for (i in 1:n) {
    S_  <- i;            V_  <- n + i;      F_  <- 2*n + i
    N_  <- 3*n + i;      VF_ <- 4*n + i;    NV_ <- 5*n + i
    NF_ <- 6*n + i
    S7_ <- 7*n + i;      V7_ <- 8*n + i;    F7_ <- 9*n + i
    N7_ <- 10*n + i;     VF7_<- 11*n + i;   NV7_<- 12*n + i
    NF7_<- 13*n + i
    m <- function(vac_idx, unv_idx, wf) {
      dm <- s[vac_idx] * wf
      s[vac_idx] <<- s[vac_idx] - dm
      s[unv_idx] <<- s[unv_idx] + dm
    }
    m(S7_, S_,   wfV); m(V7_, V_,   wfV); m(N7_, N_,   wfV); m(NV7_, NV_, wfV)
    m(F7_, F_,   wfF); m(VF7_, VF_, wfF); m(NF7_, NF_, wfF)
  }
  
  #mortality + aging + births, Euler substeps (dt = 1/n_substeps yr)
  for (sub in seq_len(n_substeps)) {
    pop_a <- numeric(n)
    for (k in 0:13) pop_a <- pop_a + s[(k * n + 1):((k + 1) * n)]
    total  <- sum(pop_a)
    births <- birth_rate_y * total
    
    #uniform mortality + aging across all 14 carriage families per age
    for (k in 0:13) {
      idx  <- (k * n + 1):((k + 1) * n)
      orig <- s[idx]
      new  <- orig * (1 - (aging_rate + mu_age) * dt)
      new[2:n] <- new[2:n] + orig[1:(n - 1)] * aging_rate[1:(n - 1)] * dt
      s[idx] <- new
    }
    
    #births enter S<1y / S_7<1y only
    s[1]         <- s[1]         + (1 - vacc_cov) * births * dt
    s[7 * n + 1] <- s[7 * n + 1] +        vacc_cov  * births * dt
  }
  
  #reset incidence accumulators
  for (k in 14:19) s[(k * n + 1):((k + 1) * n)] <- 0
  s
}


symmetrise_R <- function(c_raw_daily, pop_vec, n_age) {
  Cs <- matrix(0, n_age, n_age)
  for (i in 1:n_age) {
    for (j in 1:n_age) {
      Cs[i, j] <- (c_raw_daily[i, j]
                   + c_raw_daily[j, i] * pop_vec[j] / pop_vec[i]) / 2
    }
  }
  Cs * 365.25
}


#extract posterior draws
draws_df <- posterior::as_draws_df(fit$draws())
n_draws  <- nrow(draws_df)
N_PROJ_DRAWS <- min(N_PROJ_DRAWS, n_draws)
set.seed(1)
draw_idx <- sample(seq_len(n_draws), N_PROJ_DRAWS)

#terminal state (end of 2009, 80 columns per draw)
yes_cols <- sprintf("year_end_state[%d,%d]", n_year_stan, seq_len(n_states))
missing_cols <- setdiff(yes_cols, names(draws_df))
if (length(missing_cols)) {
  stop("Cannot find year_end_state in fit. Missing: ",
       paste(head(missing_cols, 3), collapse = ", "), " ...")
}

yes_terminal_mat <- as.matrix(draws_df[, yes_cols])

delta_V_draws <- draws_df$delta_V
delta_F_draws <- draws_df$delta_F
omega_V_draws <- draws_df$omega_V
omega_F_draws <- draws_df$omega_F
q_age_mat  <- as.matrix(draws_df[, sprintf("q_age[%d]", seq_len(n_age))])
rr_V_draws <- draws_df$rr_V
rr_F_draws <- draws_df$rr_F
rr_N_draws <- draws_df$rr_N

#ccr_1999 [age, serotype]. outer() fills a first then k, so column j corresponds to (a = ((j-1) %% n_age) + 1, k = ((j-1) %/% n_age) + 1).
ccr_cols <- as.vector(outer(seq_len(n_age), 1:3, function(a, k) sprintf("ccr_1999[%d,%d]", a, k)))
ccr_mat <- as.matrix(draws_df[, ccr_cols])

#per-draw forward simulation
#boundary rates: Stan applies age_vaccinate_demog once more at the end of its last year (using birth_rate_year[n_year] and vacc_cov_year[n_year]) to
#move the state from end-of-2009 to start-of-2010, but does not export the result. We reproduce that step below with the same rates.
br_terminal <- as.numeric(birth_rate[as.character(max(all_years))])
br_proj_prev <- as.numeric(birth_rate[as.character(PROJ_YEARS - 1)])
vacc_cov_terminal <- 0.95

pred_proj <- 
  array(NA_real_, 
        dim = c(N_PROJ_DRAWS, length(PROJ_YEARS), n_age, 3),
        dimnames = list(NULL, as.character(PROJ_YEARS), age_labels, serotypes)
        )

message(sprintf("projecting %d draws forward through %d ...(takes a few minutes)", N_PROJ_DRAWS, PROJ_END_YEAR))
pb <- txtProgressBar(min = 0, max = N_PROJ_DRAWS, style = 3)

for (d_i in seq_len(N_PROJ_DRAWS)) {
  d       <- draw_idx[d_i]
  q_age_d <- q_age_mat[d, ]
  rr_V_d  <- rr_V_draws[d]; rr_F_d <- rr_F_draws[d]; rr_N_d <- rr_N_draws[d]
  dV_d    <- delta_V_draws[d]; dF_d <- delta_F_draws[d]
  oV_d    <- omega_V_draws[d]; oF_d <- omega_F_draws[d]
  ccr_d   <- matrix(ccr_mat[d, ], nrow = n_age, ncol = 3) #[age, serotype]
  
  #move from end-of-2009 to start-of-2010
  cur <- 
    age_vaccinate_demog_R(yes_terminal_mat[d, ], 
                          n_age, aging_rate, 
                          mu_age, 
                          vacc_cov_terminal, 
                          oV_d, 
                          oF_d, 
                          br_terminal)
  
  for (t_i in seq_along(PROJ_YEARS)) {
    
    #FOI + contact symmetrisation from the compartment-summed pop
    pop_now <- numeric(n_age)
    for (k in 0:13) pop_now <- pop_now + cur[(k * n_age + 1):((k + 1) * n_age)]
    c_sym <- symmetrise_R(contact_raw, pop_now, n_age)
    
    parms <- list(n_age = n_age, q_age = q_age_d,
                  gamma_V = gamma_V, gamma_F = gamma_F, gamma_N = gamma_N,
                  rr_V = rr_V_d, rr_F = rr_F_d, rr_N = rr_N_d,
                  delta_V = dV_d, delta_F = dF_d,
                  c_sym = c_sym, pop_j = pop_now)
    
    sol <- deSolve::ode(y = cur, times = c(0, 1),
                        func = pneumo_rhs_R, parms = parms,
                        method = "lsoda", rtol = ODE_RTOL, atol = ODE_ATOL)
    endv <- as.numeric(sol[2, -1])
    
    #per-age incidence -> pred_ipd using this draw's CCR_1999
    for (a in seq_len(n_age)) {
      inc_Va <- endv[14 * n_age + a] + endv[17 * n_age + a]
      inc_Fa <- endv[15 * n_age + a] + endv[18 * n_age + a]
      inc_Na <- endv[16 * n_age + a] + endv[19 * n_age + a]
      # pred_proj[d_i, t_i, a, 1] <- inc_Va * ccr_d[a, 1]
      # pred_proj[d_i, t_i, a, 2] <- inc_Fa * ccr_d[a, 2]
      # pred_proj[d_i, t_i, a, 3] <- inc_Na * ccr_d[a, 3]
      
      pred_proj[d_i, t_i, a, 1] <- rpois(1, pmax(inc_Va * ccr_d[a, 1], 0) + 1e-6)
      pred_proj[d_i, t_i, a, 2] <- rpois(1, pmax(inc_Fa * ccr_d[a, 2], 0) + 1e-6)
      pred_proj[d_i, t_i, a, 3] <- rpois(1, pmax(inc_Na * ccr_d[a, 3], 0) + 1e-6)
    }
    
    #boundary step for the next year
    cur <- age_vaccinate_demog_R(endv, n_age, aging_rate, mu_age,
                                 VACC_COV_PROJ[t_i], oV_d, oF_d,
                                 br_proj_prev[t_i])
  }
  setTxtProgressBar(pb, d_i)
}
close(pb)


#summarise the results
proj_long <- 
  as.data.frame.table(pred_proj, responseName = "pred_ipd") %>%
  as_tibble() %>%
  setNames(c("draw_i", "year", "age", "serotype", "pred_ipd")) %>%
  mutate(year = as.integer(as.character(year)))
proj_long

proj_summary <- 
  proj_long %>%
  dplyr::group_by(year, age, serotype) %>%
  dplyr::summarise(lo = quantile(pred_ipd, 0.025, na.rm = TRUE), 
            med = quantile(pred_ipd, 0.5, na.rm = TRUE), 
            hi = quantile(pred_ipd, 0.975, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(source = "R projection")

proj_summary

#combine with Stan's pred_ipd_rep for 1999..n_year_stan
stan_pred <- 
  posterior::as_draws_df(fit$draws("pred_ipd_rep", format = "draws_array")) %>%
  posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(t)) %>%
  dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(t >= 1, t <= n_year_stan, a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  transmute(year = all_years[t], age = age_labels[a], serotype = serotypes[k], lo = `2.5%`, med = `50%`, hi = `97.5%`, source = "Stan (1999-2009)")

all_pred <- 
  bind_rows(stan_pred, proj_summary) %>%
  mutate(age = factor(age, levels = age_labels),
         serotype = factor(serotype, levels = c("F", "V", "N")))

obs_df <- 
  as.data.frame.table(obs_ipd_array(), responseName = "obs") %>% #obs_ipd_array_ext() for 1999-2019
  dplyr::rename(year = Var1, age = Var2, serotype = Var3) %>%
  mutate(year = as.integer(as.character(year)),
         age  = factor(age, levels = age_labels),
         serotype = factor(serotype, levels = c("F", "V", "N")))

oos_shade <- data.frame(xmin = max(all_years) + 0.5, xmax = PROJ_END_YEAR + 0.5)

p_proj <- 
  ggplot(data = all_pred, aes(x = year)) + 
  geom_rect(data = oos_shade, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey85", alpha = 0.5) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = serotype), alpha = 0.25) +
  geom_line(aes(y = med, colour = serotype), size =1) +
  geom_point(data = obs_df, aes(y = obs), size = 1.5, stroke = 1.2, shape = 1) +
  geom_vline(xintercept = max(all_years) + 0.5, linetype = "dashed", colour = "grey30") +
  facet_wrap(factor(serotype, levels=c('F','V','N')) ~ factor(age, levels=c('<1y', '1-4y', '5-17y', '18+y')), scales = 'free') +
  scale_x_continuous(breaks = seq(1999, PROJ_END_YEAR, by = 4)) +
  labs(title = sprintf("", PROJ_END_YEAR), y = "IPD cases", x = "Year") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

p_proj
ggsave(here::here("output", "ipd_model","forward_projection_plot.png"), p_proj, width = 16, height = 12, dpi = 150)
write.csv(proj_summary, here::here("output", "ipd_model", "projection_posterior.csv"), row.names = FALSE)
