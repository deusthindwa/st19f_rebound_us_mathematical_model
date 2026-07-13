# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(tidytable)
  library(patchwork)
})

#import datasets
source(here::here('script', "ipd_model", "01_data_prep.R"))
fit <- readRDS(here::here('results', "ipd_model", "ipd_fit.rds"))
inputs <- readRDS(here::here('results', "ipd_model", "ipd_fit_inputs.rds"))
all_years <- inputs$all_years # 1999..2009

# out-of-sample projections, falls back to n_year if not present.
n_fit <- if (!is.null(inputs$n_fit)) inputs$n_fit else length(all_years)
fit_end_year <- if (!is.null(inputs$fit_end_year)) inputs$fit_end_year else all_years[n_fit]
oos_years <- if (n_fit < length(all_years)) all_years[(n_fit + 1):length(all_years)] else integer(0)

#independent log-uniform prior bounds for the uncertainty-propagation
#parameters (set in 04_fit_model.R from carriage-fit marginal quantiles).
log_q_lower <- inputs$log_q_lower
log_q_upper <- inputs$log_q_upper
log_rr_V_lb <- inputs$log_rr_V_lb; log_rr_V_ub <- inputs$log_rr_V_ub
log_rr_F_lb <- inputs$log_rr_F_lb; log_rr_F_ub <- inputs$log_rr_F_ub
log_rr_N_lb <- inputs$log_rr_N_lb; log_rr_N_ub <- inputs$log_rr_N_ub
n_age       <- fixed_data$n_age
age_labels  <- fixed_data$age_labels
serotypes   <- c("V", "F", "N")

fit_years <- all_years  #all years in the likelihood
pcv_years <- all_years[all_years >= 2000] #2000..2009 (post-PCV)
foi_years <- 2000:2007 #FOI plot focus

draws_array <- fit$draws()
draws_df    <- posterior::as_draws_df(draws_array)

#posterior summary (natural-scale estimated parameters)
print(fit$summary(c("delta_V", "delta_F", "omega_V", "omega_F", "duration_V", "duration_F")))

#posterior summary (log-scale sampled parameters)
print(fit$summary(c("log_delta_V", "log_delta_F", "log_omega_V", "log_omega_F", "log_q_age", "log_rr_V", "log_rr_F", "log_rr_N")))

#posterior summary (uncertainty-propagation parameters)
print(fit$summary(c("q_age", "rr_V", "rr_F", "rr_N")))

#trace plots
p_trace <- 
  bayesplot::mcmc_trace(draws_array, pars = c("delta_V", "delta_F", "omega_V", "omega_F", "duration_V", "duration_F")) + 
  ggtitle("Trace plots (natural scale): chain mixing across all chains")

p_trace
ggsave(here::here('output', "ipd_model", "trace_plot.png"), p_trace, width = 9, height = 7, dpi = 150)


p_trace_log <- 
  bayesplot::mcmc_trace(draws_array, pars = c("log_delta_V", "log_delta_F", "log_omega_V", "log_omega_F")) + 
  ggtitle("Trace plots (log scale, as sampled): chain mixing across all chains")

p_trace_log
ggsave(here::here('output', "ipd_model", "traceLog_plot.png"), p_trace_log, width = 9, height = 6, dpi = 150)

#pairs plot
p_pairs <- 
  bayesplot::mcmc_pairs(draws_array, pars = c("delta_V", "delta_F", "omega_V", "omega_F"), off_diag_args = list(size = 0.3, alpha = 0.4))

p_pairs
ggsave(here::here('output', "ipd_model", "pairs_plot.png"), p_pairs, width = 8, height = 8, dpi = 150)

#joint posterior density pairs plot
#use unconstrained log-scale draws for the pairs plot
#space HMC actually explored, where correlations are most informative
param_names <- c("delta_V", "delta_F", "omega_V", "omega_F")
post_pairs <-
  as_draws_df(fit$draws(param_names)) %>%
  dplyr::select(all_of(param_names))

lower_density <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density_2d_filled(alpha = 0.85, bins = 8, show.legend = FALSE) +
    geom_point(alpha = 0.05, size = 0.4) +
    scale_fill_viridis_d(option = "mako", direction = -1) +
    theme_minimal(base_size = 8)
}

diag_density <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density(fill = "steelblue", alpha = 0.5, colour = "steelblue4") +
    theme_minimal(base_size = 8)
}

upper_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r <- cor(x, y)
  ggplot(data, mapping) +
    geom_density_2d(colour = "steelblue4", linewidth = 0.4) +
    annotate("text", x = mean(range(x)), y = mean(range(y)), label = sprintf("r = %.2f", r), size = 3, colour = "firebrick", fontface = "bold") +
    theme_minimal(base_size = 8)
}

pairs_density_plot <-
  ggpairs(post_pairs,
          lower = list(continuous = lower_density),
          diag  = list(continuous = diag_density),
          upper = list(continuous = upper_cor),
          axisLabels = "show",
          title = "") + #Joint posterior density of the 7 free log-scale parameters
  theme(strip.text = element_text(size = 7))

print(pairs_density_plot)
ggsave(here::here("output", "ipd_model", "pairs_density_advanced_plot.png"), pairs_density_plot, width = 14, height = 14, dpi = 330)

#prior vs posterior density (estimated parameters)
#priors on the natural-scale variables (sampled on the log scale, Jacobian applied in the Stan model block)
#delta_V, delta_F ~ Beta(13, 7)  equivalent to mean 0.70, 95 % CrI ~ 0.488-0.874
#omega_V, omega_F ~ Gamma(shape = 2, rate = 8) equivalent to  mean 0.25, 95 % CrI ~ 0.030-0.697
prior_df <- 
  bind_rows(
    tibble(par = "delta_V", x = seq(0.01, 0.99, length.out = 400), density = dbeta(x, 13, 7)),
    tibble(par = "delta_F", x = seq(0.01, 0.99, length.out = 400), density = dbeta(x, 13, 7)),
    tibble(par = "omega_V", x = seq(0.01, 0.99, length.out = 400), density = dgamma(x, shape = 11, rate = 4)), #2,8
    tibble(par = "omega_F", x = seq(0.01, 0.99, length.out = 400), density = dgamma(x, shape = 11, rate = 4)))

post_long <- 
  draws_df %>%
  as_tibble() %>%
  dplyr::select(delta_V, delta_F, omega_V, omega_F) %>%
  tidyr::pivot_longer(everything(), names_to = "par", values_to = "x")

p_dens <-
  ggplot() +
  geom_density(data = post_long, aes(x = x, fill = "Posterior"), alpha = 0.6, colour = NA) +
  geom_line(data = prior_df, aes(x = x, y = density, colour = "Prior"), size = 0.8) +
  facet_wrap(~par, scales = "free", ncol = 2) +
  labs(x = NULL, y = "density", title = "Prior vs posterior (natural scale)", subtitle = "Sampled on the log scale; delta ~ Beta(14, 6), omega ~ Gamma(2, 8); bounds = prior 95% CrI") +
  scale_fill_manual(name = NULL, values = c(Posterior = "#2c7bb6")) +
  scale_colour_manual(name = NULL, values = c(Prior = "black")) +
  theme_bw()

p_dens
ggsave(here::here('output', "ipd_model", "prior_posterior_plot.png"), p_dens, width = 8, height = 6, dpi = 150)

#prior vs posterior for the uncertainty-propagation parameters
#q_age and rr_V/F/N have independent log-uniform priors on the natural
#scale, with bounds equal to the marginal quantiles (default 0.005/0.995)
#of the carriage-fit posterior. Density of a log-uniform on (lo, hi) is
#d(x) = 1 / (x * (log(hi) - log(lo))) for lo <= x <= hi, else 0.
phi_names <- c(paste0("q_age[", seq_len(n_age), "]"), "rr_V", "rr_F", "rr_N")
nat_lo_vec <- c(exp(log_q_lower), exp(log_rr_V_lb), exp(log_rr_F_lb), exp(log_rr_N_lb))
nat_hi_vec <- c(exp(log_q_upper), exp(log_rr_V_ub), exp(log_rr_F_ub), exp(log_rr_N_ub))

nuis_long <- 
  draws_df %>% as_tibble() %>%
  dplyr::select(starts_with("q_age"), rr_V, rr_F, rr_N) %>%
  pivot_longer(everything(), names_to = "par", values_to = "x")

log_unif_density <- function(x, lo, hi) {
  ifelse(x >= lo & x <= hi, 1 / (x * (log(hi) - log(lo))), 0)
}

prior_nuis <- bind_rows(lapply(seq_along(phi_names), function(j) {
  lo <- nat_lo_vec[j]; hi <- nat_hi_vec[j]
  pad <- 0.05 * (hi - lo)
  xs  <- seq(max(lo - pad, .Machine$double.eps), hi + pad, length.out = 400)
  tibble(par = phi_names[j],
         x = xs,
         density = log_unif_density(xs, lo, hi))
}))

p_dens_nuis <- 
  ggplot() +
  geom_density(data = nuis_long, aes(x = x, fill = "Posterior"), alpha = 0.6, colour = NA) +
  geom_line(data = prior_nuis, aes(x = x, y = density, colour = "Prior"), size = 0.7) +
  facet_wrap(~par, scales = "free", ncol = 3) +
  labs(x = NULL, y = "density", title = "Uncertainty-propagation parameters: prior vs posterior", subtitle = "Prior = independent log-uniform between marginal quantiles of carriage_fit.rds") +
  scale_fill_manual(name = NULL, values = c(Posterior = "#d7191c")) +
  scale_colour_manual(name = NULL, values = c(Prior = "black")) +
  theme_bw()

p_dens_nuis
ggsave(here::here('output', "ipd_model", "nuisance_prior_posterior_plot.png"), p_dens_nuis, width = 9, height = 6, dpi = 150)

#FOI by year/age/serotype for 2000-2007
foi_draws <- fit$draws("foi_year", format = "draws_array")
foi_summary <- 
  posterior::summarise_draws(posterior::as_draws_df(foi_draws), ~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(t), !is.na(a), !is.na(k)) %>%
  dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(t >= 1, t <= length(all_years), a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  dplyr::mutate(year = all_years[t], age = age_labels[a], serotype = serotypes[k])

foi_plot_data <- foi_summary

p_foi <- 
  ggplot(foi_plot_data, aes(x = year, y = `50%`, colour = serotype)) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = serotype), alpha = 0.25, colour = NA) +
  facet_wrap(~ age) +
  scale_x_continuous(breaks = foi_years) +
  labs(title = "Posterior FOI by year, age and serotype (2000-2007)", y = "FOI (per year)") +
  theme_bw()

p_foi
ggsave(here::here('output', "ipd_model", "foi_year_plot.png"), p_foi, width = 9, height = 6, dpi = 150)
write.csv(foi_summary, here::here('output', "ipd_model", "foi_year_posterior_estimates.csv"), row.names = FALSE)

#CCRs in 1999 (posterior summary from generated quantities)
inc = 100000
ccr_draws <- fit$draws("ccr_1999", format = "draws_array")
ccr_summary <- 
  posterior::as_draws_df(ccr_draws) %>%
  posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(a), !is.na(k)) %>%
  dplyr::mutate(a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  dplyr::mutate(age = age_labels[a], serotype = serotypes[k], age = factor(age, levels = age_labels))
ccr_summary

p_ccr1999 <- 
  ggplot(ccr_summary,aes(x = age, y = `50%`*inc, ymin = `2.5%`*inc, ymax = `97.5%`*inc, colour = serotype)) +
  geom_pointrange(position = position_dodge(width = 0.4), shape = 4, stroke = 1.2, size = 1) +
  theme_bw(base_size = 14, base_family = "Lato") + 
  facet_wrap(.~age, scales = 'free', ncol = 4) +
  labs(title = "", y = "Case Carrier Ratios (CCRs) in 1999\nper 100,000, with 95% CrI", x ="") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 14))

p_ccr1999
ggsave(here::here('output', "ipd_model", "ccr_1999_posterior_plot.png"), p_ccr1999, width = 8, height = 5, dpi = 150)
write.csv(ccr_summary, here::here('output', "ipd_model", "ccr_1999_posterior_estimates.csv"), row.names = FALSE)

#implied 2009 CCRs (using posterior median carriage incidence of 2009 and IPD cases of 2010)
inc_draws <- fit$draws("inc_year", format = "draws_array")
inc_summary <- 
  posterior::as_draws_df(inc_draws) %>%
  posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(t), !is.na(a), !is.na(k)) %>%
  dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(t >= 1, t <= length(all_years), a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  dplyr::mutate(year = all_years[t], age = age_labels[a], serotype = serotypes[k]) %>%
  dplyr::filter(year == 2009)
inc_summary

ccr_2009_df <- 
  inc_summary %>%
  rowwise() %>%
  mutate(obs = obs_ipd_array_ext()["2010", age, serotype], 
         `2.5%` = obs/`2.5%`,
         `50%`  = obs/`50%`,
         `97.5%` = obs/`97.5%`) %>% 
  ungroup()

p_ccr2009 <- 
  ggplot(ccr_2009_df, aes(x = age, y = log(`50%`*inc), ymin = log(`2.5%`*inc), ymax = log(`97.5%`*inc), colour = serotype)) +
  geom_pointrange(position = position_dodge(width = 0.4), shape = 4, stroke = 1.2, size = 1) +
  theme_bw(base_size = 14, base_family = "Lato") + 
  facet_wrap(.~age, scales = 'free', ncol=4) +
  labs(title = "", y = "Case Carrier Ratios (CCRs) in 2009\nper 100,000, with 95% CrI", x ="") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 14))

p_ccr2009
print(ccr_2009_df %>% select(age, serotype, obs, `2.5%`, `50%`, `97.5%`))
write.csv(ccr_2009_df %>% select(age, serotype, obs, `2.5%`, `50%`, `97.5%`), here::here('output', "ipd_model", "ccr_2009.csv"), row.names = FALSE)

p_ccr1999 / p_ccr2009

#Observed vs predicted IPD (1999..2009)
pred_summary <- 
  posterior::as_draws_df(fit$draws("pred_ipd_rep", format = "draws_array")) %>%
  posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(t), !is.na(a), !is.na(k)) %>%
  dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(t >= 1, t <= length(all_years), a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  dplyr::mutate(year = all_years[t], age = age_labels[a], serotype = serotypes[k])
  
obs_df <- 
  as.data.frame.table(fixed_data$obs_ipd_array, responseName = "obs") %>%
  dplyr::rename(year = Var1, age = Var2, serotype = Var3) %>%
  mutate(year = as.integer(as.character(year))) %>%
  filter(year %in% all_years)

obs_pred <- left_join(pred_summary, obs_df, by = c("year", "age", "serotype"))

ci_df <- tidytable::map_dfr(obs_pred$obs, function(x) {
  ci <- poisson.test(x, conf.level = 0.95)$conf.int
  tibble(lower_95CI = ci[1], upper_95CI = ci[2])
})

obs_pred <- 
  bind_cols(obs_pred, ci_df) %>%
  dplyr::mutate(serotype = factor(serotype, levels = c("F", "V", "N")))

#shade out-of-sample projection window so the fit boundary is visually obvious.
oos_shade <- if (length(oos_years) > 0) {
  data.frame(xmin = fit_end_year + 0.5, xmax = max(all_years) + 0.5)
} else NULL

p_obs_pred <-
ggplot(obs_pred, aes(x = year)) +
  
  {if (!is.null(oos_shade))
  geom_rect(data = oos_shade, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey85", alpha = 0.5)} +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = serotype), alpha = 0.25) +
  geom_line(aes(y = `50%`, colour = serotype)) +
  geom_point(aes(y = obs), size = 1.2) +
  geom_errorbar(aes(x = year, y = obs, ymin = lower_95CI, ymax = upper_95CI), color = 'black', width = 0, size = 0.5) +
  geom_vline(xintercept = 1999.5, linetype = "dotted", colour = "grey50") +
  {if (length(oos_years) > 0)
  geom_vline(xintercept = fit_end_year + 0.5, linetype = "dashed", colour = "grey30")} +
  facet_wrap(factor(serotype, levels=c('F','V','N')) ~ factor(age, levels=c('<1y', '1-4y', '5-17y', '18+y')), scales = 'free') +
  labs(title = "", x = 'Year', y = "IPD cases") +
  theme_bw(base_size = 14, base_family = "Lato")
  
p_obs_pred
ggsave(here::here('output', "ipd_model", "obs_pred_plot.png"), p_obs_pred, width = 10, height = 8, dpi = 150)


#population trajectory & population-weighted FOI
pop_long <- 
  as_tibble(fixed_data$pop_matrix, rownames = "year") %>%
  mutate(year = as.integer(year)) %>%
  pivot_longer(-year, names_to = "age", values_to = "pop") %>%
  mutate(age = factor(age, levels = age_labels))

p_pop <- 
  ggplot(pop_long, aes(x = year, y = pop, colour = age)) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = foi_years) +
  labs(title = "US population by age, 2000-2007", y = "Persons") +
  theme_bw()

p_pop
ggsave(here::here('output', "ipd_model", "population.png"), p_pop, width = 8, height = 5, dpi = 150)

foi_overall <- 
  foi_plot_data %>%
  dplyr::left_join(pop_long, by = c("year", "age")) %>%
  dplyr::group_by(year, serotype) %>%
  dplyr::summarise(foi_avg = weighted.mean(`50%`, w = pop, na.rm = TRUE), .groups = "drop")

p_foi_overall <- 
  ggplot(foi_overall, aes(x = year, y = foi_avg, colour = serotype)) +
  geom_line(size = 0.8) +
  scale_x_continuous(breaks = foi_years) +
  labs(title = "Population-weighted FOI by serotype (2000-2007)", y = "FOI (per year, weighted)") +
  theme_bw()

p_foi_overall
ggsave(here::here('output', "ipd_model", "foi_overall_plot.png"), p_foi_overall, width = 8, height = 5, dpi = 150)


#out-of-sample validation for years > fit_end_year
#these years are NOT in the Poisson likelihood, pred_ipd_rep is a genuine forward projection. 
#if the observed dots sit inside the 95% predictive interval, calibration generalizes and if not the model has drift the likelihood didn't see
ppc_pred_int <- 
  posterior::as_draws_df(fit$draws("pred_ipd_rep", format = "draws_array")) %>%
  posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
  tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
  dplyr::filter(!is.na(t), !is.na(a), !is.na(k)) %>%
  dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k)) %>%
  dplyr::filter(t >= 1, t <= length(all_years), a >= 1, a <= n_age, k >= 1, k <= length(serotypes)) %>%
  dplyr::mutate(year = all_years[t], age = age_labels[a], serotype = serotypes[k])

  ppc_df <- ppc_pred_int %>%
    dplyr::filter(year %in% oos_years) %>%
    dplyr::left_join(obs_df, by = c("year", "age", "serotype"))
  
  ci_df <- map_dfr(ppc_df$obs, function(x) {
    ci <- poisson.test(x, conf.level = 0.95)$conf.int
    tibble(lower_95CI = ci[1], upper_95CI = ci[2])
  })
  
  ppc_df <- 
    dplyr::bind_cols(ppc_df, ci_df) %>% 
    dplyr::mutate(serotype = factor(serotype, levels = c("F", "V", "N")))
  
  p_ppc <- 
    ggplot(ppc_df, aes(x =factor(year))) +
    geom_pointrange(aes(y = `50%`, ymin = `2.5%`, ymax = `97.5%`, colour = serotype), position = position_dodge2(width = 1.5), size = 0.5) +
    geom_point(aes(y = obs, shape = "Observed"), colour = "black", size = 2) +
    geom_errorbar(aes(x = as.factor(year), y = obs, ymin = lower_95CI, ymax = upper_95CI), color = 'black', width = 0, size = 0.5) +
    facet_wrap(factor(serotype, levels=c('F','V','N')) ~ factor(age, levels=c('<1y', '1-4y', '5-17y', '18+y')), scales = 'free') +
    labs(title = "", y = "IPD cases", x = 'Year') +
    theme_bw(base_size = 14, base_family = "Lato") +
    theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
    theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
    theme(legend.position = 'right') +
    theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14))
  
  p_ppc
  ggsave(here::here('output', "ipd_model", "ppc_out_of_sample_plot.png"), p_ppc, width = 16, height = 12, dpi = 150)


#carriage prevalence over time (V, F, N), 1999-2009
#carriage in <1y develops during the year via the FOI, so the meaningful measurement is at year END.
#we compute year-end carriage from year_end_state in Stan
#unvaccinated S, V, F, N, VF, NV, NF; and vaccinated S_7, V_7, F_7, N_7, VF_7, NV_7, NF_7
n_states_per_age <- 20 # 14 carriage + 6 incidence
fam_V   <- c(1, 4, 5, 8, 11, 12)
fam_F   <- c(2, 4, 6, 9, 11, 13)
fam_N   <- c(3, 5, 6, 10, 12, 13)
fam_all <- 0:13 
n_year_loc <- length(all_years)

#prev_year_end, already exported by Stan model
has_prev_end <- tryCatch({
  d <- fit$draws("prev_year_end", format = "draws_array")
  !is.null(d) && length(d) > 0
}, error = function(e) FALSE)


if (has_prev_end) {
  prev_summary <- 
    posterior::as_draws_df(fit$draws("prev_year_end", format = "draws_array")) %>%
    posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
    tidyr::extract(variable, into = c("var", "t", "a", "k"), regex = "^([A-Za-z_0-9]+)\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$", remove = TRUE) %>%
    dplyr::filter(!is.na(t), !is.na(a), !is.na(k)) %>%
    dplyr::mutate(t = as.integer(t), a = as.integer(a), k = as.integer(k))
} else {
  message("prev_year_end not in this fit; reconstructing from year_end_state.")
  yes_mat <- posterior::as_draws_matrix(fit$draws("year_end_state", format = "draws_array"))
  n_draws  <- nrow(yes_mat)
  col_name <- function(t, idx) sprintf("year_end_state[%d,%d]", t, idx)

  prev_end <- array(NA_real_,
                    dim = c(n_draws, n_year_loc, n_age, 3))
  for (t in seq_len(n_year_loc)) {
    for (a in seq_len(n_age)) {
      pop_cols  <- vapply(fam_all,
                          function(k) col_name(t, k * n_age + a), "")
      pop_age_a <- rowSums(yes_mat[, pop_cols, drop = FALSE])
      for (k_idx in 1:3) {
        fam       <- list(fam_V, fam_F, fam_N)[[k_idx]]
        carr_cols <- vapply(fam,
                            function(k) col_name(t, k * n_age + a), "")
        carr      <- rowSums(yes_mat[, carr_cols, drop = FALSE])
        prev_end[, t, a, k_idx] <- carr / pmax(pop_age_a, 1e-9)
      }
    }
  }

  prev_summary <- expand.grid(t = seq_len(n_year_loc),
                              a = seq_len(n_age),
                              k = 1:3) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      `2.5%`  = quantile(prev_end[, t, a, k], 0.025, na.rm = TRUE),
      `50%`   = quantile(prev_end[, t, a, k], 0.5,   na.rm = TRUE),
      `97.5%` = quantile(prev_end[, t, a, k], 0.975, na.rm = TRUE)
    ) %>% ungroup()
}

prev_summary <- 
  prev_summary %>%
  mutate(year = all_years[t], age = age_labels[a], serotype = serotypes[k]) %>%
  mutate(age = factor(age, levels = age_labels)) %>%
  mutate(serotype = factor(serotype, levels = c("F", "V", "N")))

p_prev <- 
  ggplot(prev_summary, aes(x = year, y = `50%`, colour = serotype, fill = serotype)) +
  geom_point(size = 1.8, shape = 1, stroke = 1.5) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.25, colour = NA) +
  geom_vline(xintercept = 1999.5, linetype = "dotted", colour = "grey50") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~ age, scales = "free_y") +
  scale_x_continuous(breaks = seq(1999, 2009, by = 1)) +
  labs(title = "", y = "Carriage prevalence", x = "Year") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14))

p_prev
prev_summary
ggsave(here::here('output', "ipd_model", "carriage_prevalence_plot.png"), p_prev, width = 10, height = 6, dpi = 150)
write.csv(prev_summary %>% select(year, age, serotype, `2.5%`, `50%`, `97.5%`), here::here('output', "ipd_model", "carriage_prevalence_posterior.csv", row.names = FALSE))

#summary tables
param_summary <- 
  fit$summary(c("delta_V", "delta_F", "omega_V", "omega_F","duration_V", "duration_F", "q_age", "rr_V", "rr_F", "rr_N"),
  mean = mean, 
  median = median,
  q2.5 = ~ quantile(.x, 0.025), 
  q97.5 = ~ quantile(.x, 0.975),
  rhat = posterior::rhat, 
  ess_bulk = posterior::ess_bulk)

param_summary
write.csv(param_summary, here::here('output', "ipd_model", "parameter_summary.csv"), row.names = FALSE)
