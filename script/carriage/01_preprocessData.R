# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#estimates of carriage from English raw data from UKHSA
#English data processing
df_EW <-
  rio::import(here::here('data', 'carriage', 'eng_carr_studies.xlsx')) %>%
  dplyr::filter(yearc == '2001/2') %>%
  dplyr::mutate(STg = if_else(is.na(ST), NA_character_,
                              if_else(ST == '19F', '19F',
                                      if_else(grepl("\\b(4|6A|6B|9V|14|18C|23F)\\b", ST) == TRUE, 'VT', 'NVT'))),
                Age_groups = factor(if_else(Age <1, "<1y",
                                     if_else(Age >=1 & Age <5, "1-4y",
                                             if_else(Age >=5 & Age <18, "5-17y",
                                                     if_else(Age >=18, "18+y", NA_character_)))), levels = c("<1y","1-4y","5-17y","18+y"))) %>%
  dplyr::group_by(STg, Age_groups) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STg) %>%
  tidyr::pivot_wider(values_from = c(n), names_from = STg) %>%
  dplyr::ungroup()

df_EW <-
  df_EW %>%
  dplyr::rowwise() %>%
  dplyr::mutate(N.carr = sum(`19F`, VT, NVT),
                N.total = sum(`19F`, VT, NVT, `NA`),
                Uncol = `NA`/N.total,
                Car.prev = N.carr/N.total,
                VT.prev = VT/N.total,
                F.prev = `19F`/N.total,
                NVT.prev = NVT/N.total) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Population = c(629200, 3009200, 10131900, 39274600))

df_EWo <-
  df_EW %>%
  dplyr::select(Age_groups, VT.prev, F.prev, NVT.prev, N.carr, N.total, Population) %>%
  tidyr::pivot_longer(cols = c(VT.prev, F.prev, NVT.prev)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Age_groups, name, N.carr, N.total, value) %>%
  dplyr::mutate(n = N.total*value) %>%
  dplyr::mutate(Low = 1/(1 + (N.total - n + 1)/(n * qf(0.05, 2 * n, 2 * (N.total - n + 1)))),
                High = 1/(1 + (N.total - n)/((n + 1) * qf(1 - 0.05, 2 * (n + 1), 2 * (N.total - n)))),
                carrtype = "Observed",
                stg = if_else(name == 'VT.prev', 'V',
                              if_else(name == 'NVT.prev', 'N', 'F'))) %>%
  dplyr::rename('agegp'="Age_groups", 'Est'='value') %>%
  dplyr::select(agegp, Est, Low, High, stg, carrtype)

#====================================================================

# #total contact matrix 
# age_groups_5yr <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
# 
# #weighted mean daily in-person contacts (Epistorm-Mix, 2024)
# M <- matrix(c(
#   1.86068,0.70922,0.27294,0.19284,0.38647,0.52152,0.74264,0.65879,0.50745,0.27078,0.13497,0.11624,0.20563,0.22665,0.07370,0.03396,
#   0.21983,4.18650,0.88506,0.37158,0.19047,0.35853,0.66981,0.76950,0.98274,0.32328,0.37390,0.25236,0.15199,0.14732,0.11631,0.08913,
#   0.16316,0.24448,3.97903,0.65597,0.29650,0.11857,0.39405,0.54558,0.78884,0.61383,0.39497,0.11900,0.14986,0.08006,0.07873,0.06910,
#   0.28803,0.66286,0.52931,3.81800,0.63310,0.36918,0.47484,0.51248,0.78651,0.77182,0.60021,0.29266,0.15755,0.13396,0.11359,0.16432,
#   0.24193,0.21859,0.51182,0.80917,1.22247,0.63708,0.48638,0.65223,0.69792,0.49391,0.50489,0.48721,0.35503,0.11931,0.11213,0.19379,
#   0.39399,0.39625,0.16909,0.32176,0.45123,0.95707,0.81098,0.34533,0.35295,0.24431,0.49553,0.34436,0.34734,0.13904,0.08809,0.17117,
#   0.41114,0.40073,0.61230,0.47785,0.54435,0.66417,1.14187,0.95482,0.60610,0.38127,0.50529,0.46032,0.48470,0.27719,0.14004,0.15737,
#   0.60228,0.57665,1.03357,0.39680,0.35656,0.45246,0.74558,0.98167,0.64406,0.59712,0.39807,0.30514,0.35053,0.27646,0.25654,0.17135,
#   0.28746,0.38422,0.48765,0.77529,0.46967,0.51555,0.51592,0.60443,0.93232,0.77716,0.35524,0.46906,0.39576,0.22935,0.29411,0.21664,
#   0.09829,0.36638,0.57683,0.97212,0.37505,0.35279,0.48063,0.60579,0.76518,0.89284,0.77723,0.52809,0.40137,0.21828,0.13325,0.26808,
#   0.19709,0.26375,0.19000,0.39801,0.57346,0.61348,0.86635,0.61171,0.49417,0.62418,0.86114,0.70865,0.51590,0.17011,0.16663,0.18595,
#   0.32640,0.25551,0.48136,0.58423,0.71343,0.70075,0.71431,0.68880,0.70071,0.62075,0.76345,1.05835,0.53343,0.25145,0.17170,0.36944,
#   0.21101,0.17483,0.12477,0.12894,0.27466,0.39818,0.60666,0.33268,0.43297,0.32044,0.38106,0.37968,0.68928,0.33411,0.20188,0.25076,
#   0.17283,0.15241,0.20782,0.17603,0.23082,0.30466,0.39359,0.43447,0.30474,0.27823,0.34361,0.26130,0.46794,0.76394,0.39551,0.36527,
#   0.03191,0.05615,0.12077,0.08700,0.09648,0.15435,0.25678,0.18681,0.32028,0.33532,0.41090,0.24662,0.50743,0.47868,0.65688,0.54496,
#   0.00000,0.02814,0.06460,0.06333,0.11065,0.26669,0.30602,0.17459,0.33583,0.30125,0.33912,0.28585,0.26830,0.41705,0.38301,1.15239
# ), nrow = 16, byrow = TRUE,
# dimnames = list(respondent = age_groups_5yr, contact = age_groups_5yr))
# 
# #population sizes (single-year)
# single_yr <- c(
#   `0`=3480117, `1`=3532512, `2`=3672703, `3`=3797741, `4`=3917162,
#   `5`=4001330, `6`=3975522, `7`=4018926, `8`=4059908, `9`=4074737,
#   `10`=4251732,`11`=4268172,`12`=4421435,`13`=4388774,`14`=4297717,
#   `15`=4322527,`16`=4343161,`17`=4281824,`18`=4417941,`19`=4670623)
# 
# #population for each 5-yr group used in the matrix
# N <- 
#   list(
#     "lt1"  = sum(single_yr[c("0")]),
#     "1to4" = sum(single_yr[c("1","2","3","4")]),
#     "0to4" = sum(single_yr[as.character(0:4)]),
#     "5to9" = sum(single_yr[as.character(5:9)]),
#     "10to14" = sum(single_yr[as.character(10:14)]),
#     "15to17" = sum(single_yr[c("15","16","17")]),
#     "18to19" = sum(single_yr[c("18","19")]),
#     "15to19" = sum(single_yr[as.character(15:19)]),
#     # Approximate 2024 Census estimates for ages 20+
#     "20to24" = 21200000, "25to29" = 22800000, "30to34" = 23500000,
#     "35to39" = 22300000, "40to44" = 21000000, "45to49" = 20100000,
#     "50to54" = 21500000, "55to59" = 22700000, "60to64" = 21100000,
#     "65to69" = 18700000, "70to74" = 16500000, "75plus" = 24300000
#   )
# 
# N_18plus <- 
#   N[["18to19"]] + N[["20to24"]] + N[["25to29"]] + N[["30to34"]] + N[["35to39"]] + N[["40to44"]] + N[["45to49"]] + 
#   N[["50to54"]] + N[["55to59"]] + N[["60to64"]] + N[["65to69"]] + N[["70to74"]] +  N[["75plus"]]
# 
# N_5to17 <- 
#   N[["5to9"]] + N[["10to14"]] + N[["15to17"]]
# 
# #fractional weights for splitting mixed 5-yr groups
# w_lt1_in04    <- N[["lt1"]]    / N[["0to4"]]
# w_1to4_in04   <- N[["1to4"]]  /  N[["0to4"]]
# w_15to17_in19 <- N[["15to17"]] / N[["15to19"]]
# w_18to19_in19 <- N[["18to19"]] / N[["15to19"]]
# 
# #aggregate ROWS to broad respondent age groups
# # M_broad[I, ] = Σ_{i∈I} (N_i / N_I) × M[i, ]
# 
# #<1y and 1-4y, both mapped from 0-4 row (finest resolution available)
# r_lt1   <- M["0-4", ]
# r_1to4  <- M["0-4", ]
# 
# #5-17y: pop-weighted average of 5-9, 10-14, and the 15-17 portion of 15-19
# r_5to17 <- (N[["5to9"]]   * M["5-9",   ] + N[["10to14"]] * M["10-14", ] + N[["15to17"]] * M["15-19", ]) / N_5to17
# 
# #18+y: pop-weighted average of 18-19 portion of 15-19, plus 20-24 through 75+
# age20plus_groups <- c("20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
# N_20plus_vec <- c(N[["20to24"]], N[["25to29"]], N[["30to34"]], N[["35to39"]],
#                   N[["40to44"]], N[["45to49"]], N[["50to54"]], N[["55to59"]],
#                   N[["60to64"]], N[["65to69"]], N[["70to74"]], N[["75plus"]])
# 
# r_18plus <- (N[["18to19"]] * M["15-19", ] + colSums(N_20plus_vec * M[age20plus_groups, ])) / N_18plus
# 
# #aggregate columns to broad contact age groups
# #for contact group J: sum M[i, j] for j∈J
# #(gives total contacts respondent i has with all individuals in group J)
# col_agg <- function(row_vec) {
#   c(`<1y`   = w_lt1_in04   * row_vec["0-4"],
#     `1-4y`  = w_1to4_in04  * row_vec["0-4"],
#     `5-17y` = row_vec["5-9"] + row_vec["10-14"] + w_15to17_in19 * row_vec["15-19"],
#     `18+y`  = w_18to19_in19 * row_vec["15-19"] + sum(row_vec[age20plus_groups]))
# }
# 
# M_broad <- rbind(
#   `<1y`   = col_agg(r_lt1),
#   `1-4y`  = col_agg(r_1to4),
#   `5-17y` = col_agg(r_5to17),
#   `18+y`  = col_agg(r_18plus)
# )
# colnames(M_broad) <- c("<1y","1-4y","5-17y","18+y")
# 
# #ymmetrise
# #C_ij = (M_ij + M_ji × N_j / N_i) / 2
# pop_broad <- c(`<1y`   = N[["lt1"]],
#                `1-4y`  = N[["1to4"]],
#                `5-17y` = N_5to17,
#                `18+y`  = N_18plus)
# 
# broad_labels <- c("<1y","1-4y","5-17y","18+y")
# 
# C_sym <- matrix(NA, nrow = 4, ncol = 4,
#                 dimnames = list(respondent = broad_labels,
#                                 contact    = broad_labels))
# 
# for (i in broad_labels) {
#   for (j in broad_labels) {
#     C_sym[i, j] <- (M_broad[i, j] + M_broad[j, i] * pop_broad[j] / pop_broad[i]) / 2
#   }
# }
# 
# contact_EW <- M_broad

#contact_EW <- read.csv(curl("https://raw.githubusercontent.com/StefanFlasche/Pneumo_Trans_Inf/refs/heads/master/Data/AvNumberOfContacts_EW_new.csv"))
#contacts from Flasche et al. compressed into 4 age groups
#daily average number of physical contacts between any two specific individuals in age i and age j
contact_EW <-
  as.matrix(rbind(c(0.40, 0.06, 0.02, 0.04),
                  c(0.80, 1.88, 0.33, 0.43),
                  c(0.53, 0.74, 2.24, 0.52),
                  c(2.37, 1.87, 1.38, 1.42))
            )
  