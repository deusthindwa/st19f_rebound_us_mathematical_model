# by Deus Thindwa, Dan Weinberger, Ginny Pitzer
# age-structured mathematical model for pneumococcal transmission
# 05/23/2026

#====================================================================

suppressPackageStartupMessages({
  library(tibble)
})

#age groups
age_labels <- c("<1y", "1-4y", "5-17y", "18+y")
n_age <- length(age_labels)

#years spent in each age group = aging rate (alpha_i = 1 / years)
years_in_age <- c(1, 4, 13, 59)
aging_rate <- 1 / years_in_age  # per year

#carriage susceptibility (per infectious contact)
#bounds for uncertainty propagation, the model draws q from Uniform(q_lower, q_upper) PER MCMC iteration (it is not fit to data).
q_lower <- c(`<1y`   = 0.0125, `1-4y`  = 0.0150, `5-17y` = 0.00864, `18+y`  = 0.00517)
q_upper <- c(`<1y`   = 0.0174, `1-4y`  = 0.0185, `5-17y` = 0.0108, `18+y`  = 0.00635)
q_mean <- (q_lower + q_upper)/2 #used only for the lightweight R seed

#carriage clearance rates per year
carriage_dur_days <- c(71, 34, 18, 17)
gamma_age <- 365.25/carriage_dur_days

#carriage of V, F and N all clear with the same age-specific rate.
gamma_V <- gamma_age
gamma_F <- gamma_age
gamma_N <- gamma_age

#relative risks of co-colonisation
#bounds for uncertainty propagation (Stan samples each rr on log scale). The mean is used only by the R seed burn-in to produce a near-steady-state seed state.
rr_V_lower <- 0.0712; rr_V_upper <- 0.198
rr_F_lower <- 0.352;  rr_F_upper <- 0.952
rr_N_lower <- 0.178;  rr_N_upper <- 0.524
rr_V <- (rr_V_lower + rr_V_upper)/2
rr_F <- (rr_F_lower + rr_F_upper)/2
rr_N <- (rr_N_lower + rr_N_upper)/2

#PCV7 coverage in infants (<1y), value not used in Stan, just for lightweight R seed
vacc_cov_infants <- 0.86 

#unsymmetrised contact matrix
#daily average contacts: rows = participant age, columns = contact age.
contact_raw <- matrix(c(
  0.352, 1.509, 0.865, 3.819,
  0.352, 1.509, 0.865, 3.819,
  0.055, 0.236, 4.582, 4.833,
  0.050, 0.215, 0.999, 5.702), nrow = 4, byrow = TRUE,
dimnames = list(age_labels, age_labels))

#symmetrise the contact matrix for a given year's age-group population sizes.
#returns daily contacts, multiply by 365.25 for annual.
symmetrise_contacts <- function(C_raw, pop_vec) {
  n  <- nrow(C_raw)
  Cs <- matrix(0, n, n, dimnames = dimnames(C_raw))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      Cs[i, j] <- (C_raw[i, j] + C_raw[j, i] * pop_vec[j] / pop_vec[i]) / 2
    }
  }
  Cs
}

#observed US population 1999-2019
us_pop <- tibble::tribble(
  ~year,    ~`<1y`,    ~`1-4y`,   ~`5-17y`,  ~`18+y`,
  1999,   258310,  1033240,  3497789, 14271440,
  2000,   331435,  1325739,  4031650, 16470672,
  2001,   335992,  1343966,  4065801, 16694393,
  2002,   366761,  1467046,  4532325, 18609391,
  2003,   370365,  1481461,  4540649, 18801641,
  2004,   400515,  1602061,  4908700, 20392829,
  2005,   414130,  1656520,  4916630, 20650525,
  2006,   416014,  1664058,  4956165, 20984777,
  2007,   420870,  1683480,  4983002, 21282436,
  2008,   422993,  1691973,  5011344, 21565745,
  2009,   436450,  1745798,  4968817, 21395015,
  2010,   427053,  1708210,  5114878, 21983399,
  2011,   427391,  1709562,  5559073, 22269108,
  2012,   426224,  1704898,  5569599, 22551687,
  2013,   424061,  1696244,  5580304, 22820875,
  2014,   422180,  1688721,  5589490, 23096981,
  2015,   420923,  1683693,  5597785, 23706554,
  2016,   419468,  1677874,  5604720, 23952680,
  2017,   418554,  1674214,  5614354, 26009514,
  2018,   417798,  1671192,  5616201, 26254610,
  2019,   416608,  1666432,  5612473, 26499965
)
pop_matrix <- as.matrix(us_pop[, -1])  # rows = year, cols = age group
rownames(pop_matrix) <- us_pop$year
pop_matrix2 <- pop_matrix

#observed IPD case counts 1999..2019
obs_ipd <- tibble::tribble(
  ~year, ~V_1, ~V_2, ~V_3, ~V_4, ~F_1, ~F_2, ~F_3, ~F_4, ~N_1, ~N_2, ~N_3, ~N_4,
  1999,   597,  163,   61, 1510,   82,   24,    6,  107,   91,   28,   49, 1066,
  2000,   546,  155,   79, 1555,   61,   18,    6,  112,  103,   31,   37, 1069,
  2001,   208,  156,   64, 1433,   21,   37,    8,   96,  129,   45,   63, 1225,
  2002,    60,   74,   61, 1108,   19,   12,    9,   74,  145,   64,   59, 1465,
  2003,    55,   27,   54,  822,   17,    6,    5,   71,  187,   48,   60, 1744,
  2004,    22,   13,   27,  651,   10,   10,    7,   54,  203,   52,   67, 1847,
  2005,    12,    8,   21,  482,    8,    7,    1,   77,  236,   57,   95, 2322,
  2006,     7,    5,   15,  367,    5,    3,    3,   46,  240,   63,  124, 2524,
  2007,     4,    4,    6,  272,    4,    0,    3,   30,  277,   73,  121, 2668,
  2008,     4,    6,    4,  226,    1,    2,    1,   23,  257,   58,  117, 2930,
  2009,     2,    2,    7,  159,    1,    2,    1,   26,  277,   80,  147, 2958,
  2010,   141,   96,   65, 1340,    1,    1,    1,   26,   89,   54,   40, 1599, #EXCEPTIONAL: the number for F_3 cases in 2010 was 0 but replaced it with 1 so that CCRs are well defined for 2010-2019
  2011,    36,   39,   49, 1069,    3,    1,    2,   20,  111,   52,   45, 1781,
  2012,    15,   15,   33,  709,    1,    2,    5,   32,   92,   55,   57, 1757,
  2013,     9,   18,   21,  619,    2,    0,    5,   37,   93,   58,   54, 2058,
  2014,    12,   13,   10,  479,    4,    4,    1,   41,   77,   57,   59, 1792,
  2015,    13,    8,   12,  497,    7,    3,    4,   62,   88,   42,   45, 1832,
  2016,    11,    8,   21,  540,    5,    6,    6,   74,   74,   41,   43, 1910,
  2017,    12,   11,   19,  624,    5,    4,   11,   72,   55,   41,   54, 2028,
  2018,     8,    9,   20,  637,    5,    3,   11,   72,   70,   38,   52, 1996,
  2019,     9,   11,   14,  597,    9,    2,    3,   88,   68,   38,   49, 1963
)

#return array [n_year, n_age, n_serotype], serotype order = c("V","F","N")
obs_ipd_array <- function() {
  yrs  <- obs_ipd$year
  arr  <- array(NA_real_,
                dim = c(length(yrs), n_age, 3),
                dimnames = list(yrs, age_labels, c("V", "F", "N")))
  cols <- c("V_1","V_2","V_3","V_4","F_1","F_2","F_3","F_4","N_1","N_2","N_3","N_4")
  m    <- as.matrix(obs_ipd[, cols])
  for (k in seq_len(3)) {
    arr[, , k] <- m[, ((k - 1) * 4 + 1):(k * 4)]
  }
  arr
}

#bundle everything
fixed_data <- list(
  age_labels       = age_labels,
  n_age            = n_age,
  aging_rate       = aging_rate,
  q_lower          = q_lower,    q_upper    = q_upper,    q_mean = q_mean,
  rr_V_lower       = rr_V_lower, rr_V_upper = rr_V_upper, rr_V_mean = rr_V,
  rr_F_lower       = rr_F_lower, rr_F_upper = rr_F_upper, rr_F_mean = rr_F,
  rr_N_lower       = rr_N_lower, rr_N_upper = rr_N_upper, rr_N_mean = rr_N,
  rr_V             = rr_V, rr_F = rr_F, rr_N = rr_N,
  gamma_V          = gamma_V,
  gamma_F          = gamma_F,
  gamma_N          = gamma_N,
  vacc_cov_infants = vacc_cov_infants,
  contact_raw      = contact_raw,
  pop_matrix       = pop_matrix,
  obs_ipd          = obs_ipd,
  obs_ipd_array    = obs_ipd_array()
)

#PCV7 / PCV13 parameter ranges (each scenario passes its own bounds for the PCV7- and PCV13-specific values).
#common across all 3 scenarios
delta_V_7_lower <- 0.68;  delta_V_7_upper   <- 0.75
delta_V_13_lower <- 0.68;  delta_V_13_upper  <- 0.75
omega_V_7_lower <- 0.10;  omega_V_7_upper   <- 0.20
omega_V_13_lower <- 0.10;  omega_V_13_upper  <- 0.20
omega_F_7_lower <- 0.10;  omega_F_7_upper   <- 0.20

#delta_F_7 (PCV7 efficacy against F)
delta_F_7_lower   <- 0.65
delta_F_7_upper   <- 0.70

#delta_F_13 (PCV13 efficacy against F):
delta_F_13_lower_prop <- 0.65
delta_F_13_upper_prop <- 0.70

#omega_F_13 (PCV13 waning against F):
omega_F_13_lower_prop <- 0.10
omega_F_13_upper_prop <- 0.20

#re-bundle fixed_data with everything
fixed_data$pop_matrix <- pop_matrix
fixed_data$pop_matrix2 <- pop_matrix2
fixed_data$obs_ipd <- obs_ipd
fixed_data$obs_ipd_array <- obs_ipd_array()

fixed_data$delta_V_7_lower <- delta_V_7_lower;  fixed_data$delta_V_7_upper <- delta_V_7_upper
fixed_data$delta_V_13_lower <- delta_V_13_lower; fixed_data$delta_V_13_upper <- delta_V_13_upper
fixed_data$omega_V_7_lower <- omega_V_7_lower;  fixed_data$omega_V_7_upper <- omega_V_7_upper
fixed_data$omega_V_13_lower <- omega_V_13_lower; fixed_data$omega_V_13_upper <- omega_V_13_upper
fixed_data$omega_F_7_lower <- omega_F_7_lower;  fixed_data$omega_F_7_upper <- omega_F_7_upper
fixed_data$delta_F_7_lower <- delta_F_7_lower;  fixed_data$delta_F_7_upper <- delta_F_7_upper
fixed_data$delta_F_13_lower_prop <- delta_F_13_lower_prop
fixed_data$delta_F_13_upper_prop <- delta_F_13_upper_prop
fixed_data$omega_F_13_lower_prop <- omega_F_13_lower_prop
fixed_data$omega_F_13_upper_prop <- omega_F_13_upper_prop
