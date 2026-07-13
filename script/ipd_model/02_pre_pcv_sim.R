# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

suppressPackageStartupMessages({
  library(deSolve)
})

#run data preparation script
source(here::here("script", "ipd_model", "01_data_prep.R"))
#source(here::here("01_data_prep.R"))

#compartment index
# State vector layout (length 28 + 12 = 40):
#   1- 4 : S    (age 1..4)
#   5- 8 : V
#   9-12 : F
#  13-16 : N
#  17-20 : VF
#  21-24 : NV
#  25-28 : NF
#  29-32 : cumulative incidence of V
#  33-36 : cumulative incidence of F
#  37-40 : cumulative incidence of N
N_COMP <- 7
N_INC <- 3
COMP_NMS <- c("S", "V", "F", "N", "VF", "NV", "NF")
INC_NMS <- c("incV", "incF", "incN")

idx_comp <- function(comp) {
  k <- match(comp, COMP_NMS)
  ((k - 1) * fixed_data$n_age + 1):(k * fixed_data$n_age)
}

idx_inc  <- function(comp) {
  k <- match(comp, INC_NMS)
  (N_COMP * fixed_data$n_age + (k - 1) * fixed_data$n_age + 1):
    (N_COMP * fixed_data$n_age + k * fixed_data$n_age)
}

state_init <- function(target_pop) {
  s0 <- rep(0, N_COMP * fixed_data$n_age + N_INC * fixed_data$n_age)
  s0[idx_comp("S")]  <- target_pop * 0.94
  s0[idx_comp("V")]  <- target_pop * 0.03
  s0[idx_comp("F")]  <- target_pop * 0.005
  s0[idx_comp("N")]  <- target_pop * 0.025
  s0
}

pre_pcv_ode <- function(t, y, parms) {
  S  <- y[idx_comp("S")];  V  <- y[idx_comp("V")];  Fc <- y[idx_comp("F")]
  Nc <- y[idx_comp("N")];  VF <- y[idx_comp("VF")]; NV <- y[idx_comp("NV")]
  NF <- y[idx_comp("NF")]

  carr_V <- V + VF + NV
  carr_F <- Fc + VF + NF
  carr_N <- Nc + NV + NF

  pop_j   <- parms$pop_j
  c_sym_a <- parms$c_sym_annual
  q       <- parms$q_age

  lamV <- q * as.vector(c_sym_a %*% (carr_V / pop_j))
  lamF <- q * as.vector(c_sym_a %*% (carr_F / pop_j))
  lamN <- q * as.vector(c_sym_a %*% (carr_N / pop_j))

  gV <- parms$gamma_V; gF <- parms$gamma_F; gN <- parms$gamma_N
  rV <- parms$rr_V;   rF <- parms$rr_F;   rN <- parms$rr_N

  dS  <- -S * (lamV + lamF + lamN) + gV * V + gF * Fc + gN * Nc
  dV  <-  S * lamV - V * (gV + rV * (lamF + lamN)) + gF * VF + gN * NV
  dF  <-  S * lamF - Fc * (gF + rF * (lamV + lamN)) + gV * VF + gN * NF
  dN  <-  S * lamN - Nc * (gN + rN * (lamV + lamF)) + gV * NV + gF * NF
  dVF <-  V * rV * lamF + Fc * rF * lamV - VF * (gV + gF)
  dNV <-  V * rV * lamN + Nc * rN * lamV - NV * (gV + gN)
  dNF <-  Nc * rN * lamF + Fc * rF * lamN - NF * (gN + gF)

  d_incV <- S * lamV + rF * Fc * lamV + rN * Nc * lamV
  d_incF <- S * lamF + rV * V  * lamF + rN * Nc * lamF
  d_incN <- S * lamN + rV * V  * lamN + rF * Fc * lamN

  out <- numeric(length(y))
  out[idx_comp("S")]  <- dS
  out[idx_comp("V")]  <- dV
  out[idx_comp("F")]  <- dF
  out[idx_comp("N")]  <- dN
  out[idx_comp("VF")] <- dVF
  out[idx_comp("NV")] <- dNV
  out[idx_comp("NF")] <- dNF
  out[idx_inc("incV")] <- d_incV
  out[idx_inc("incF")] <- d_incF
  out[idx_inc("incN")] <- d_incN
  list(out)
}

year_boundary <- function(y, target_pop) {
  n_age <- fixed_data$n_age
  out   <- y
  for (comp in COMP_NMS) {
    idx <- idx_comp(comp)
    x   <- y[idx]
    new <- x
    for (i in seq_len(n_age)) {
      out_rate <- fixed_data$aging_rate[i]
      new[i]   <- x[i] * (1 - out_rate)
      if (i > 1) new[i] <- new[i] + x[i - 1] * fixed_data$aging_rate[i - 1]
    }
    out[idx] <- new
  }
  out[idx_comp("S")[1]] <- out[idx_comp("S")[1]] + target_pop[1]
  totals <- sapply(seq_len(n_age), function(i) {
    sum(sapply(COMP_NMS, function(c) out[idx_comp(c)[i]]))
  })
  for (comp in COMP_NMS) {
    idx <- idx_comp(comp)
    for (i in seq_len(n_age)) {
      if (totals[i] > 0) out[idx[i]] <- out[idx[i]] * (target_pop[i] / totals[i])
    }
  }
  out[idx_inc("incV")] <- 0
  out[idx_inc("incF")] <- 0
  out[idx_inc("incN")] <- 0
  out
}

#default (fixed_data$pop_matrix["1999", ]) but (04_fit_model.R) passes calibrated 1999 pop (demog$init_1999_fitted)
simulate_pre_pcv_seed <- function(n_burn_years = 40, verbose = FALSE, pop_99 = NULL) {
  n_age   <- fixed_data$n_age
  if (is.null(pop_99)) pop_99 <- fixed_data$pop_matrix["1999", ]
  c_sym_d <- symmetrise_contacts(fixed_data$contact_raw, pop_99)
  c_sym_a <- c_sym_d * 365.25

  parms <- list(
    pop_j        = pop_99,
    c_sym_annual = c_sym_a,
    q_age        = fixed_data$q_mean,
    gamma_V      = fixed_data$gamma_V,
    gamma_F      = fixed_data$gamma_F,
    gamma_N      = fixed_data$gamma_N,
    rr_V         = fixed_data$rr_V,
    rr_F         = fixed_data$rr_F,
    rr_N         = fixed_data$rr_N)
  
  state <- state_init(pop_99)
  
  for (yi in seq_len(n_burn_years)) {
    out <- ode(y = state, times = c(0, 1), func = pre_pcv_ode,
               parms = parms, method = "lsoda", rtol = 1e-8, atol = 1e-10)
    end <- as.numeric(out[2, -1])
    state <- year_boundary(end, pop_99)
    if (verbose && (yi %% 10 == 0 || yi == n_burn_years)) {
      message(sprintf("  R seed burn-in year %d", yi))
    }
  }
  state
}

if (sys.nframe() == 0) {
  seed_state <- simulate_pre_pcv_seed(verbose = TRUE)
  totals <- sapply(seq_len(fixed_data$n_age), function(i) {
    sum(sapply(COMP_NMS, function(c) seed_state[idx_comp(c)[i]]))
  })
  message("\n seed totals by age (should match 1999 table):")
  for (i in seq_len(fixed_data$n_age)) {
    message(sprintf("  %5s state=%.0f  table=%.0f", fixed_data$age_labels[i], totals[i], fixed_data$pop_matrix["1999", i]))
  }
  saveRDS(seed_state, file = here::here("results", "ipd_model", "seed_state.rds"))
  #saveRDS(seed_state, file = here::here("seed_state.rds"))
}
