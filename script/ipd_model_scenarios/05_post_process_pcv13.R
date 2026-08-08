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
  library(deSolve)
  library(rio)
  library(scales)
  has_loo <- requireNamespace("loo", quietly = TRUE)
})


#per-scenario diagnostics + cross-scenario model comparison + two new plots

#load data preparations
source(here::here("script", "ipd_model_scenarios","01_data_prep.R"))

#fit <- readRDS(here::here('results', 'ipd_model_scenarios', "ipd_fit_scenario1.rds"))
#inputs_fit <- readRDS(here::here('results', 'ipd_model_scenarios', "ipd_fit_inputs_scenario1.rds"))

#input parameters
n_age <- fixed_data$n_age
age_labels <- fixed_data$age_labels
serotypes <- c("V", "F", "N")

#tidy quantile function
.tidy_quantiles <- function(draws_array, .var) {
  posterior::as_draws_df(draws_array) %>%
    posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE), .args = list(names = c("lo", "med", "hi"))) %>%
    tidyr::separate(variable, into = c("var", "t", "a", "k"), sep = "[\\[,\\]]", extra = "drop", fill = "right")
}

#====================================================================

#per-scenario processing begins here
process_scenario <- function(scenario_id, estimated_par_name, prior_overlay) { 
  
  fit_file <- here::here('results', 'ipd_model_scenarios', sprintf("ipd_fit_scenario%d.rds", scenario_id))
  inputs_file <- here::here('results', 'ipd_model_scenarios', sprintf("ipd_fit_inputs_scenario%d.rds", scenario_id))
  if (!file.exists(fit_file) || !file.exists(inputs_file)) { message("skipping scenario ", scenario_id, " (no fit file found)") 
    return(invisible(NULL)) 
    }

  fit <- readRDS(fit_file)
  inputs <- readRDS(inputs_file)
  all_years <- inputs$all_years
  prefix <- sprintf("scenario%d", scenario_id)

  #trace plots
  pars_for_trace <- intersect(c(estimated_par_name, paste0("log_", estimated_par_name)), dimnames(fit$draws())[[3]])
  
  if (length(pars_for_trace) > 0) {
    p_trace <- bayesplot::mcmc_trace(fit$draws(), pars = pars_for_trace) +
      ggtitle(sprintf("Trace -- scenario %d (%s)", scenario_id, estimated_par_name))
    #ggsave(here::here('output', 'ipd_model_scenarios', sprintf("%s_trace_plot.png", prefix)), p_trace, width = 9, height = 5, dpi = 150)
  }
  
  #pairs plot
  post_pairs <- as_draws_df(fit$draws(estimated_par_name)) %>% dplyr::select(all_of(estimated_par_name))
  rio::export(post_pairs, here::here('results', 'ipd_model_scenarios', sprintf("%s_pairs_plot.csv", prefix)))
 
  #prior vs posterior on the estimated parameter
  post_x <- as.numeric(posterior::as_draws_df(fit$draws(estimated_par_name))[[estimated_par_name]])
  
  if (length(post_x) > 0 && !is.null(prior_overlay)) {
    p_dens <- 
      ggplot() +
      geom_density(aes(x = post_x, fill = "Posterior"), alpha = 0.6, colour = NA) +
      geom_line(aes(x = prior_overlay$grid, y = prior_overlay$density, colour = "Prior"), size = 0.8) +
      scale_fill_manual(name = NULL, values = c(Posterior = "#2c7bb6")) +
      scale_colour_manual(name = NULL, values = c(Prior = "black")) +
      labs(x = estimated_par_name, y = "density", title = sprintf("Scenario %d: prior vs posterior (%s)", scenario_id, estimated_par_name)) +
      theme_bw()
    #ggsave(here::here('output', 'ipd_model_scenarios', sprintf("%s_prior_post_plot.png", prefix)), p_dens, width = 7, height = 4, dpi = 150)
  }


  #CCR_2010 posterior
  # ccr_summary <- 
  #   .tidy_quantiles(fit$draws("ccr_2010", format = "draws_array")) %>%
  #   dplyr::mutate(t = as.integer(t), 
  #                 a = as.integer(a), 
  #                 age = age_labels[t], 
  #                 serotype = serotypes[a]) %>%
  #   dplyr::select(-var, -k)
  # 
  # write.csv(ccr_summary, here::here('output', 'ipd_model_scenarios', sprintf("%s_ccr_2010_posterior.csv", prefix)), row.names = FALSE)

  
  #observed vs predicted IPD 1999..2019
  pred_summary <- 
    .tidy_quantiles(fit$draws("pred_ipd_rep", format = "draws_array")) %>%
    dplyr::mutate(t = as.integer(t), 
                  a = as.integer(a), 
                  k = as.integer(k),
                  year = all_years[t], 
                  age = age_labels[a], 
                  serotype = serotypes[k])
  
  obs_df <- 
    as.data.frame.table(fixed_data$obs_ipd_array, responseName = "obs") %>%
    dplyr::rename(year = Var1, age = Var2, serotype = Var3) %>%
    dplyr::mutate(year = as.integer(as.character(year))) %>%
    dplyr::filter(year %in% all_years)
  
  obs_pred <- 
    dplyr::left_join(pred_summary, obs_df, by = c("year", "age", "serotype")) %>%
    dplyr::mutate(period = case_when(year <= 2009 ~ "pre-PCV13 (informational)", TRUE ~ "PCV13 fit (2010-2019)"))
  
  p_obs_pred <- 
    ggplot(obs_pred, aes(x = year)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = serotype, group = period), alpha = 0.25) +
    geom_line(aes(y = `50%`, colour = serotype, group = period, linetype = period), size = 1.5) +
    geom_point(aes(y = obs, shape = period), size = 2.5, stroke = 1.5) +
    geom_vline(xintercept = c(1999.5, 2009.5), linetype = "dotted", colour = "grey50", size = 1) +
    scale_shape_manual(values = c("pre-PCV13 (informational)" = 1, "PCV13 fit (2010-2019)"     = 16)) +
    facet_wrap(factor(serotype, levels=c('F','V','N')) ~ factor(age, levels=c('<1y', '1-4y', '5-17y', '18+y')), scales = 'free') +
    #labs(title = sprintf("Scenario %d: observed vs posterior-predicted IPD", scenario_id), subtitle = "Likelihood is over 2010-2019 only") +
    labs(x = 'Year', y = "IPD cases") +
    scale_x_continuous(limits = c(1999, 2019), breaks = seq(1999, 2019, by = 4)) +
    theme_bw(base_size = 14, base_family = "Lato") +
    theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
    theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
    theme(legend.position = 'bottom') +
    theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14))

  
  rio::export(obs_pred, here::here('results', 'ipd_model_scenarios', sprintf("%s_obs_pred_plot.csv", prefix)))
  ggsave(here::here('output', 'ipd_model_scenarios', sprintf("%s_obs_pred_plot.png", prefix)), p_obs_pred, width = 20, height = 12, dpi = 330)


  #carriage prevalence by age and serotype, 2010-2019 year_start_state is (n_year, n_states); carriers per serotype are the
  #sum of single + dual states involving that serotype, across all strata.
  yss_arr <- fit$draws("year_end_state", format = "draws_array")
  ys_df <- 
    posterior::as_draws_df(yss_arr) %>%
    posterior::summarise_draws(~ quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE)) %>%
    dplyr::rename(lo = `2.5%`, med = `50%`, hi = `97.5%`) %>%
    tidyr::extract(variable, into  = c("var", "t", "idx"), regex = "^(.+?)\\[(\\d+),(\\d+)\\]$") %>%
    dplyr::mutate(t = as.integer(t), idx = as.integer(idx))
  
  #reconstruct carriers per (year, age, serotype)
  v_fam <- c(1, 4, 5, 8, 11, 12, 15, 18, 19)   #V single + dual + PCV strata
  f_fam <- c(2, 4, 6, 9, 11, 13, 16, 18, 20)
  n_fam <- c(3, 5, 6, 10, 12, 13, 17, 19, 20)
  
  # build_prev <- function(med, fams, age_i, t) {
  #   sapply(fams, function(k) {
  #     med[(t - 1) * 96 + k * n_age + age_i + 1]  # 1-indexed in draws
  #   }) |> sum()
  # }
  
  ys_med  <- ys_df %>% dplyr::select(t, idx, med) %>% dplyr::arrange(t, idx)
  ys_medL <- ys_df %>% dplyr::select(t, idx, lo) %>% dplyr::arrange(t, idx)
  ys_medH <- ys_df %>% dplyr::select(t, idx, hi) %>% dplyr::arrange(t, idx)
  
  yss_mat  <- matrix(ys_med$med, nrow = max(ys_med$t), byrow = TRUE)
  yss_matL <- matrix(ys_medL$lo, nrow = max(ys_medL$t), byrow = TRUE)
  yss_matH <- matrix(ys_medH$hi, nrow = max(ys_medH$t), byrow = TRUE)
  
  prev_df <- expand.grid(year = all_years, age = age_labels, serotype = serotypes, stringsAsFactors = FALSE)
  prev_dfL <- expand.grid(year = all_years, age = age_labels, serotype = serotypes, stringsAsFactors = FALSE)
  prev_dfH <- expand.grid(year = all_years, age = age_labels, serotype = serotypes, stringsAsFactors = FALSE)
  
  prev_df$prev <- NA_real_
  prev_dfL$prevlo <- NA_real_
  prev_dfH$prevhi <- NA_real_
  
  for (i in seq_len(nrow(prev_df))) {
    t <- match(prev_df$year[i], all_years)
    a <- match(prev_df$age[i], age_labels)
    fams <- switch(prev_df$serotype[i], V = v_fam, F = f_fam, N = n_fam)
    carr <- sum(yss_mat[t, fams * n_age + a])
    total <- sum(yss_mat[t, (0:20) * n_age + a])
    prev_df$prev[i] <- if (total > 0) carr / total else NA_real_
  }
  
  for (i in seq_len(nrow(prev_dfL))) {
    t <- match(prev_dfL$year[i], all_years)
    a <- match(prev_dfL$age[i], age_labels)
    fams <- switch(prev_dfL$serotype[i], V = v_fam, F = f_fam, N = n_fam)
    carr <- sum(yss_matL[t, fams * n_age + a])
    total <- sum(yss_matL[t, (0:20) * n_age + a])
    prev_dfL$prevlo[i] <- if (total > 0) carr / total else NA_real_
  }

  for (i in seq_len(nrow(prev_dfH))) {
    t <- match(prev_dfH$year[i], all_years)
    a <- match(prev_dfH$age[i], age_labels)
    fams <- switch(prev_dfH$serotype[i], V = v_fam, F = f_fam, N = n_fam)
    carr <- sum(yss_matH[t, fams * n_age + a])
    total <- sum(yss_matH[t, (0:20) * n_age + a])
    prev_dfH$prevhi[i] <- if (total > 0) carr / total else NA_real_
  }


  prev_pcv13 <-
    prev_df %>% dplyr::left_join(prev_dfL) %>% dplyr::left_join(prev_dfH) %>%
    dplyr::filter(year >= 2010) %>%
    dplyr::mutate(age = factor(age, levels = age_labels))

  p_prev <- 
    ggplot(prev_pcv13, aes(x = year, y = prev, colour = serotype)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = prevlo, ymax = prevhi), alpha = 0.25, colour = NA) +
    geom_point(size = 1.4, shape = 1, stroke = 1.5) +
    facet_wrap(~ age, scales = "free_y") + 
    geom_vline(xintercept = 2009.5, linetype = "dotted", colour = "grey50") +
    labs(y = "Carriage prevalence", x = "Year") + #title = sprintf("Scenario %d: post-PCV13 carriage prevalence by age and serotype", scenario_id), 
    theme_bw(base_size = 14, base_family = "Lato") +
    theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
    theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
    theme(legend.position = 'bottom') +
    theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_continuous(limit = c(2009, 2019), breaks = seq(2010, 2019, by = 1))
  
  
  rio::export(prev_pcv13, here::here('results', 'ipd_model_scenarios', sprintf("%s_carriage_prevalence_plot.csv", prefix)))
  ggsave(here::here('output', 'ipd_model_scenarios', sprintf("%s_carriage_prevalence_plot.png", prefix)), p_prev, width = 14, height = 8, dpi = 150)
  invisible(list(fit = fit, all_years = all_years, inputs = inputs))
  
  
}


#====================================================================

#prior overlays for each scenario's estimated parameter
prior_overlays <- list(
  `1` = list(grid = seq(0.001, 0.99, length.out = 400), density = dbeta(seq(0.001, 0.99, length.out = 400), 3, 9)),
  `2` = list(grid = seq(0.001, 3.0, length.out = 400), density = dgamma(seq(0.001, 3.0, length.out = 400), shape = 3, rate = 3)),
  `3` = list(grid = seq(0.001, 0.99, length.out = 400), density = dbeta(seq(0.001, 0.99, length.out = 400), 7, 2)))

results <- list()
results[["1"]] <- process_scenario(1, "delta_F_13", prior_overlays[["1"]])
results[["2"]] <- process_scenario(2, "omega_F_13", prior_overlays[["2"]])
results[["3"]] <- process_scenario(3, "rr_V",  prior_overlays[["3"]]) #rr_N_post


#consolidated density plot
density_ds <-
  dplyr::bind_rows(
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario1_pairs_plot.csv')) %>% dplyr::mutate(scenario = "Reduced direct PCV effectiveness") %>% dplyr::rename('value' = 'delta_F_13'),
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario2_pairs_plot.csv')) %>% dplyr::mutate(scenario = "Reduced duration of PCV protection") %>% dplyr::rename('value' = 'omega_F_13'),
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario3_pairs_plot.csv')) %>% dplyr::mutate(scenario = "Increased risk of N to co-colonisation") %>% dplyr::rename('value' = 'rr_V')) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("Reduced direct PCV effectiveness", "Reduced duration of PCV protection", "Increased risk of N to co-colonisation")))

A <- 
  bayesplot::mcmc_trace(density_ds %>% dplyr::filter(scenario == 'Reduced direct PCV effectiveness') %>% dplyr::select(value)) + 
  labs(x = 'Iterations', y = "Parameter value") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(legend.position = 'none')

B <- 
  bayesplot::mcmc_trace(density_ds %>% dplyr::filter(scenario == 'Reduced duration of PCV protection') %>% dplyr::select(value)) + 
  labs(x = 'Iterations', y ="") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(legend.position = 'none')

C <- 
  bayesplot::mcmc_trace(density_ds %>% dplyr::filter(scenario == 'Increased risk of N to co-colonisation') %>% dplyr::select(value)) + 
  labs(x = 'Iterations', y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(legend.position = 'none')

D <- 
  density_ds %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(. ~ scenario, scales = 'free') +
  labs(x = 'Parameter value', y = "Density") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14))

(A|B|C) /D


#consolidated IPD plot
dplyr::bind_rows(
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario1_obs_pred_plot.csv')) %>% dplyr::mutate(scenario = "Reduced direct PCV effectiveness"),
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario2_obs_pred_plot.csv')) %>% dplyr::mutate(scenario = "Reduced duration of PCV protection"),
  rio::import(here::here('results', 'ipd_model_scenarios', 'scenario3_obs_pred_plot.csv')) %>% dplyr::mutate(scenario = "Increased risk of N to co-colonisation")) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("Reduced direct PCV effectiveness", "Reduced duration of PCV protection", "Increased risk of N to co-colonisation"))) %>%
  dplyr::filter(serotype =='F') %>%
  #dplyr::filter(scenario == "Reduced direct PCV effectiveness") %>%
  
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = serotype, group = period), alpha = 0.25) +
  geom_line(aes(y = `50%`, colour = serotype, group = period, linetype = period), size = 1.5) +
  geom_point(aes(y = obs, shape = period), size = 2.5, stroke = 1.5) +
  geom_vline(xintercept = c(1999.5, 2009.5), linetype = "dotted", colour = "grey50", size = 1) +
  scale_shape_manual(values = c("pre-PCV13 (propagation)" = 1, "PCV13 fit (2010-2019)" = 16)) +
  facet_wrap(scenario ~ factor(age, levels=c('<1y', '1-4y', '5-17y', '18+y')), scales = 'free') +
  labs(x = 'Year', y = "IPD cases") +
  scale_x_continuous(limits = c(1999, 2019), breaks = seq(1999, 2019, by = 4)) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'right') +
  theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 14))


#consolidated carriage prevalence plot
scenarios_prev <- 
  dplyr::bind_rows(
    rio::import(here::here('results', 'ipd_model_scenarios', 'scenario1_carriage_prevalence_plot.csv')) %>% dplyr::mutate(scenario = "Reduced direct PCV effectiveness"),
    rio::import(here::here('results', 'ipd_model_scenarios', 'scenario2_carriage_prevalence_plot.csv')) %>% dplyr::mutate(scenario = "Reduced duration of PCV protection"),
    rio::import(here::here('results', 'ipd_model_scenarios', 'scenario3_carriage_prevalence_plot.csv')) %>% dplyr::mutate(scenario = "Increased risk of N to co-colonisation")) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c("Reduced direct PCV effectiveness", "Reduced duration of PCV protection", "Increased risk of N to co-colonisation")))

scenarios_prev %>%
  ggplot(aes(x = year, y = prev, colour = scenario)) +
  geom_line(size = 1) +
  #geom_ribbon(aes(ymin = prevlo, ymax = prevhi, colour = serotype), alpha = 0.8) +
  geom_point(size = 2.4, shape = 1, stroke = 1.5) +
  facet_grid(factor(serotype, levels = c('F',"V",'N')) ~ factor(age, levels = c('<1y', '1-4y', '5-17y', '18+y')), scales = "free") + 
  labs(y = "Carriage prevalence", x = "Year") + 
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.text.y = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5), axis.text.y = element_text(size = 12)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(2010, 2019), breaks = seq(2010, 2019, by = 1)) +
  theme(legend.position = 'bottom')


#====================================================================


# cross-scenario WAIC + lppd
# computes WAIC and lppd for each scenario from log_lik.
# lppd (log pointwise predictive density) is the raw fit component of WAIC
# larger (less negative) lppd means better in-sample predictive fit; WAIC penalises complexity via p_waic.
#   lppd = sum_i log( (1/S) * sum_s exp(log_lik[s, i]) )
#   p_waic = sum_i var_s( log_lik[s, i] ) # effective parameters
#   WAIC    = -2 * (lppd - p_waic) # = -2 * loo::waic()$elpd_waic

#llpd
.logSumExp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
compute_lppd <- function(ll_mat) {
  S <- nrow(ll_mat)
  sum(apply(ll_mat, 2, .logSumExp)) - ncol(ll_mat) * log(S)
}

#WAIC
fits_present <- sapply(results, function(r) !is.null(r))
if (has_loo && all(fits_present)) {
  waic_list <- list()
  lppd_vec  <- numeric(3)
  for (s in c("1", "2", "3")) {
    fit    <- results[[s]]$fit
    ll_mat <- as.matrix(fit$draws("log_lik", format = "draws_matrix"))
    waic_list[[paste0("scenario", s)]] <- loo::waic(ll_mat)
    lppd_vec[as.integer(s)] <- compute_lppd(ll_mat)
  }
  comp <- loo::loo_compare(waic_list)
  print(comp)
  print(waic_list)
  
  lppd_df <- data.frame(
    scenario = paste0("scenario", 1:3),
    lppd     = lppd_vec,
    dlppd    = lppd_vec - max(lppd_vec)   #0 for the best-fitting scenario, negative for the others
  )
  
  cat("\nlog pointwise predictive density (lppd)\n")
  print(lppd_df, row.names = FALSE)
  
} else {
  message("\nskip WAIC comparison, `loo` is needed to run all 3 scenarios")
}

# #cross-scenario BIC [BIC = -2 * log(L_hat) + k * log(n)]
# # log(L_hat): posterior-maximum log-likelihood over the MCMC draws (~ maximum log-likelihood).
# # k : number of parameters informatively constrained by the fit-window data. Each scenario estimates ONE PCV13 parameter (delta_F_13, omega_F_13, or rr_N_post) from the 2010-2019 IPD counts
# #    the 8 propagated pre-PCV13 parameters and 2 "shadow" inactive PCV13 declarations are pinned by their propagation-prior bounds and are NOT informed by the current likelihood, so they do not enter k. 
# #    k = 1 in every scenario. Because k is identical across scenarios, delta-BIC comes entirely from the difference in log(L_hat).
# # n: observations in the fit window = 10 yr x 4 ages x 3 serotypes = 120. log_lik is computed for all 21 years in generated quantities but the model likelihood was only summed over t >= n_fit_start_t (t = 12..21)
# #.   we filter to match, so BIC reflects the same information the sampler used. (WAIC above does not filter, since all three scenarios share the same t < 12 log_lik contributions, the ranking is unaffected)
# # Kass & Raftery (1995) evidence against the lowest-BIC scenario:
# # 0-2 (barely worth mentioning), 2-6 (positive evidence for the lower-BIC model), 6-10 (strong), >10 (very strong)
# 
# compute_scenario_bic <- function(fit, k = 1L, fit_start_t = 12L) {
#   ll_draws  <- fit$draws("log_lik", format = "draws_matrix")
#   var_names <- colnames(ll_draws)
#   t_idx     <- as.integer(sub("^log_lik\\[(\\d+),.*$", "\\1", var_names))
#   keep      <- which(t_idx >= fit_start_t)
#   ll_fit    <- as.matrix(ll_draws[, keep, drop = FALSE])   # n_draws x n_obs
#   ll_per_draw <- rowSums(ll_fit)
#   n_obs     <- length(keep)
#   list(BIC       = -2 * max(ll_per_draw) + k * log(n_obs),
#        log_L_hat = max(ll_per_draw),
#        k         = k,
#        n_obs     = n_obs)
# }
# 
# #compute BIC for each scenario
# if (all(fits_present)) {
#   bic_rows <- lapply(c("1", "2", "3"), function(s) {
#     b <- compute_scenario_bic(results[[s]]$fit)
#     data.frame(scenario  = paste0("scenario", s),
#                log_L_hat = b$log_L_hat,
#                k         = b$k,
#                n_obs     = b$n_obs,
#                BIC       = b$BIC,
#                stringsAsFactors = FALSE)
#   })
#   bic_df <- do.call(rbind, bic_rows)
#   bic_df$dBIC <- bic_df$BIC - min(bic_df$BIC)
#   bic_df$evidence_against <- cut(
#     bic_df$dBIC,
#     breaks = c(-Inf, 2, 6, 10, Inf),
#     labels = c("(best or near-best)",
#                "positive (2-6)",
#                "strong (6-10)",
#                "very strong (>10)"),
#     right = FALSE
#   )
#   print(bic_df)
# } else {
#   message("\nskip, BIC comparison needs all 3 scenario fits")
# }
