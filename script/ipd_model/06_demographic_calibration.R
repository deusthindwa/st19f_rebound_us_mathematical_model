# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

suppressPackageStartupMessages({
  library(here)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

#run data preparation script
source(here::here("script", "ipd_model", "01_data_prep.R"))

#inputs
birth_df <- #https://wonder.cdc.gov/natality.html
  tibble::tibble(birth_rate = c(0.01195276, 0.011939196, 0.011926342, 0.011913077, 0.0119, 0.011886019, 0.011871513, 0.011856968, 0.011841408, 0.011825831, 0.011809191, 
                 0.011791992, 0.011775391, 0.01175709, 0.011738829, 0.011719383, 0.0117, 0.011679432, 0.011658562, 0.011638356, 0.011617758))

stopifnot("birth_rate" %in% names(birth_df))
stopifnot(nrow(birth_df) == 21) #1999 .. 2019
birth_years <- 1999:2019
br_yearly <- setNames(birth_df$birth_rate, as.character(birth_years))
INIT_1999 <- c(`<1y` = 258310, `1-4y` = 1033240, `5-17y` = 3497789, `18+y` = 14271440)
SIM_YEARS <- 1999:2019

#calibration target = every year of observed pop available in fixed_data$pop_matrix (intersected with SIM_YEARS). Extending the
#calibration window past 2009 -- if the pop table contains those years pulls mu away from the very-negative values needed to fit the steep
#2000-2010 growth alone, producing a flatter trajectory that also fits the 2010-2019 plateau.
CALIB_YEARS <- intersect(SIM_YEARS, as.integer(rownames(fixed_data$pop_matrix2)))
N_SUBSTEPS <- 12  #set monthly Euler substeps

#maximum allowed shift on the 1999 initial pop, per age, on the log scale. bounds on log_scale = +/- INIT_LOG_SCALE_BOUND map to a multiplicative shift of roughly +/- 35 %
#setting this to 0 reverts to the old "pin to
#INIT_1999 exactly" behaviour (and the calibration falls back to mu-only).
INIT_LOG_SCALE_BOUND <- 0.30
age_labels   <- names(INIT_1999)
aging_rate   <- setNames(c(1, 1/4, 1/13, 1/59), age_labels)
stopifnot(identical(unname(aging_rate), as.numeric(fixed_data$aging_rate)))

#one-year demographic step
#births enter age 1 at rate b * total. Each age outflows at rate alpha_i + mu_i (aging + net mortality). Age i > 1 also gains inflow from
#age i-1 at rate alpha_{i-1} * N_{i-1}. Solved via n_substeps Euler steps.
demog_one_year <- function(pop, birth_rate, mu, alpha = aging_rate, n_substeps = N_SUBSTEPS) {
  n  <- length(pop)
  dt <- 1 / n_substeps
  for (s in seq_len(n_substeps)) {
    total <- sum(pop)
    new   <- pop
    for (i in seq_len(n)) {
      out_rate <- alpha[i] + mu[i]
      new[i] <- pop[i] * (1 - out_rate * dt)
      if (i > 1) new[i] <- new[i] + pop[i - 1] * alpha[i - 1] * dt
    }
    new[1] <- new[1] + birth_rate * total * dt
    pop    <- new
  }
  pop
}

simulate_demog <- function(mu, init_pop, br_yearly, years, alpha = aging_rate) {
  n_y <- length(years)
  pop <- matrix(0, nrow = n_y, ncol = length(init_pop),
                dimnames = list(as.character(years), names(init_pop)))
  pop[1, ] <- init_pop
  for (i in seq_len(n_y - 1)) {
    yr <- years[i]
    b  <- br_yearly[as.character(yr)]
    pop[i + 1, ] <- demog_one_year(pop[i, ], b, mu, alpha)
  }
  pop
}

#calibrate mu + initial-pop shift against observed pop
#theta = (mu_<1y, mu_1-4y, mu_5-17y, mu_18+y, log_scale_<1y, log_scale_1-4y, log_scale_5-17y, log_scale_18+y)
#the simulated 1999 initial pop is INIT_1999 * exp(log_scale_a), so the optimiser can shift trajectory up or down per age (bounded). 
#freeing initial pop removes "1999 obs = sim exactly by construction" artefact, the curve is no longer pinned through one point
#usually drops mu's magnitude substantially.
obs_pop_calib <- as.matrix(fixed_data$pop_matrix2[as.character(CALIB_YEARS), ])

obj_fn <- function(theta) {
  mu        <- theta[1:n_age]
  log_scale <- theta[(n_age + 1):(2 * n_age)]
  init_now  <- INIT_1999 * exp(log_scale)
  pop_sim   <- simulate_demog(mu, init_now, br_yearly, CALIB_YEARS)
  
  #relative-error SSE so all four age groups contribute on equal footing.
  sum(((pop_sim - obs_pop_calib) / obs_pop_calib) ^ 2)
}

#optimisation
init_theta <- c(0.005, 0.001, 0.001, 0.005, rep(0, n_age))

opt <- optim(init_theta, 
             obj_fn,
             method = "L-BFGS-B",
             lower = c(rep(-Inf, n_age), rep(-INIT_LOG_SCALE_BOUND, n_age)),
             upper = c(rep( Inf, n_age), rep( INIT_LOG_SCALE_BOUND, n_age)),
             control = list(maxit = 5000, factr = 1e7))

mu_hat <- setNames(opt$par[1:n_age], age_labels)
log_scale_hat<- setNames(opt$par[(n_age + 1):(2 * n_age)], age_labels)
init_hat<- INIT_1999 * exp(log_scale_hat)

#calibration iterations and convergence
opt$counts[1]
opt$convergence

#SSE on relative-error scale
opt$value
length(obs_pop_calib)
length(CALIB_YEARS)

#calibrated NET mortality rates (per year)
print(round(mu_hat, 5))
mu_total <- sum(mu_hat * init_hat) / sum(init_hat)

#population-weighted average
mu_total

#calibrated 1999 initial population (shifted from given)
shift_pct <- (exp(log_scale_hat) - 1) * 100
for (a in seq_along(init_hat)) {
  message(sprintf("  %-6s : %10.0f  (%+5.1f %% shift from given %d)",
                  age_labels[a], init_hat[a], shift_pct[a], INIT_1999[a]))
}

if (any(mu_hat < 0)) {
  message("If one or more mu_hat values are negative, this is a net demographic rate (mortality minus net immigration)")
}

#simulate full 1999-2019 trajectory using the shifted initial pop
model_pop_matrix <- simulate_demog(mu_hat, init_hat, br_yearly, SIM_YEARS)

#persist for downstream scripts
saveRDS(list(
  mu_age           = mu_hat,
  birth_rate       = br_yearly,
  model_pop_matrix = model_pop_matrix,
  init_1999_given  = INIT_1999, #observed 1999 initial pop
  init_1999_fitted = init_hat, #after calibrated log_scale shift
  log_scale_hat    = log_scale_hat,
  age_labels       = age_labels,
  sim_years        = SIM_YEARS,
  calib_years      = CALIB_YEARS,
  calib_sse        = opt$value,
  n_substeps       = N_SUBSTEPS
), file = here::here("results", "ipd_model", "demographic_model.rds"))

#overall observed vs simulated total population plot
obs_all <- fixed_data$pop_matrix2
obs_yrs <- as.integer(rownames(obs_all))

p_overall <-
  bind_rows(tibble(year = obs_yrs, total_pop = rowSums(obs_all), type = "Observed"),
            tibble(year = SIM_YEARS, total_pop = rowSums(model_pop_matrix), type = "Model-simulated")) %>%
  mutate(type = factor(type, levels = c("Observed", "Model-simulated"))) %>%
  ggplot(aes(x = year, y = total_pop, colour = type, shape = type, group = type)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.4) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(1999, 2019, by = 2)) +
  scale_colour_manual(values = c("Observed" = "black", "Model-simulated" = "#1f77b4")) +
  scale_shape_manual(values = c("Observed" = 16, "Model-simulated" = 17)) +
  labs(title = "Overall population: observed vs model-simulated", x = "Year", y = "Total population", colour = NULL, shape = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

p_overall
ggsave(here::here("output", "ipd_model", "plot_demog_overall.png"), p_overall, width = 9, height = 5, dpi = 150)

#age-group bar plot, observed vs simulated plot
obs_long <- 
  as_tibble(obs_all, rownames = "year") %>%
  mutate(year = as.integer(year)) %>%
  pivot_longer(-year, names_to = "age", values_to = "pop") %>%
  mutate(type = "Observed")

sim_long <- 
  as_tibble(model_pop_matrix, rownames = "year") %>%
  mutate(year = as.integer(year)) %>%
  pivot_longer(-year, names_to = "age", values_to = "pop") %>%
  mutate(type = "Model")

bar_df <- bind_rows(obs_long, sim_long) %>%
  mutate(age  = factor(age, levels = age_labels), 
         type = factor(type, levels = c("Observed", "Model")))

# p_bars <- 
#   ggplot(bar_df, aes(x = year, y = pop/sum(pop), fill = type)) +
#   geom_col(position = position_dodge(width = 0.85, preserve = "single"), width = 0.75) +
#   facet_wrap(~ age, scales = "free_y") +
#   scale_y_continuous(labels = scales::comma) +
#   scale_x_continuous(breaks = seq(1999, 2019, by = 2)) +
#   scale_fill_manual(values = c("Observed" = "grey40", "Model" = "#d95f02")) +
#   labs(title = "Age-group population: observed vs model-simulated", x = "Year", y = "Proportion of population", fill = NULL) +
#   theme_bw() +
#   theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

p_bars <- 
  bar_df %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(p=pop/sum(pop)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = p, fill = type)) +
  geom_col(position = position_dodge(width = 0.85, preserve = "single"), width = 0.75) +
  facet_wrap(~ age, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(1999, 2019, by = 2)) +
  scale_fill_manual(values = c("Observed" = "grey40", "Model" = "#d95f02")) +
  labs(title = "Age-group population: observed vs model-simulated", x = "Year", y = "Proportion of population", fill = NULL) +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'bottom') +
  theme(axis.text.x = element_text(size = 14, angle = 30, hjust = 1), axis.text.y = element_text(size = 14))

p_bars
ggsave(here::here("output", "ipd_model", "plot_demog_by_age.png"), p_bars, width = 11, height = 6, dpi = 150)
