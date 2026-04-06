# by Deus Thindwa and Dan Weinberger
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#pcv7 annual cases
ipd7 <-
  rio::import(here::here("data", "usa_ipd.csv")) %>%
  dplyr::filter(st != "MISS", yearc <2010) %>%
  dplyr::mutate(agegp = if_else(agegp == 'Age <2', '0-1y', if_else(agegp == 'Age 2-4', '2-4y', if_else(agegp == 'Age 5-17', '5-17y', '18+y'))),
                yearc = as.integer(yearc),
                yearc = if_else(yearc == 1998, 1999, yearc)) %>%
  dplyr::group_by(yearc, agegp, st) %>%
  dplyr::tally(N_cases) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(agegp = factor(agegp, levels(factor(agegp))[c(1, 3, 4, 2)])) %>%
  dplyr::mutate(stg = if_else(st %in% c('4', '6A', '6B', '9V', '14', '18C', '23F'), "V", if_else(st == '19F', 'F', 'N')),
                age = if_else(agegp == '0-1y', 'gp1', if_else(agegp == '2-4y', 'gp2', if_else(agegp == '5-17y', 'gp3', 'gp4')))) %>%
  dplyr::group_by(yearc, age, stg) %>%
  dplyr::summarise(n = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = as.integer(if_else(yearc == 1999, n/2, n)))

#pcv7A annual cases
ipd7A <-
  rio::import(here::here("data", "usa_ipd.csv")) %>%
  dplyr::filter(st != "MISS", yearc >=2010, yearc <= 2019) %>%
  dplyr::mutate(agegp = if_else(agegp == 'Age <2', '0-1y', if_else(agegp == 'Age 2-4', '2-4y', if_else(agegp == 'Age 5-17', '5-17y', '18+y'))),
                yearc = as.integer(yearc),
                yearc = if_else(yearc == 1998, 1999, yearc)) %>%
  dplyr::group_by(yearc, agegp, st) %>%
  dplyr::tally(N_cases) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(agegp = factor(agegp, levels(factor(agegp))[c(1, 3, 4, 2)])) %>%
  mutate(stg = if_else(st %in% c('4', '6A', '6B', '9V', '14', '18C', '23F'), "V", if_else(st == '19F', 'F', 'N')),
         age = if_else(agegp == '0-1y', 'gp1', if_else(agegp == '2-4y', 'gp2', if_else(agegp == '5-17y', 'gp3', 'gp4')))) %>%
  dplyr::group_by(yearc, age, stg) %>%
  dplyr::summarise(n = sum(n)) %>%
  dplyr::ungroup()

#pcv13 annual cases
ipd13 <-
  rio::import(here::here("data", "usa_ipd.csv")) %>%
  dplyr::filter(st != "MISS", yearc >=2010, yearc <= 2019) %>%
  dplyr::mutate(agegp = if_else(agegp == 'Age <2', '0-1y', if_else(agegp == 'Age 2-4', '2-4y', if_else(agegp == 'Age 5-17', '5-17y', '18+y'))),
                yearc = as.integer(yearc),
                yearc = if_else(yearc == 1998, 1999, yearc)) %>%
  dplyr::group_by(yearc, agegp, st) %>%
  dplyr::tally(N_cases) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(agegp = factor(agegp, levels(factor(agegp))[c(1, 3, 4, 2)])) %>%
  mutate(stg = if_else(st %in% c('4', '6A', '6B', '9V', '14', '18C', '23F', '1', '5', '6A', '7F', '19A'), "V", if_else(st == '19F', 'F', 'N')),
         age = if_else(agegp == '0-1y', 'gp1', if_else(agegp == '2-4y', 'gp2', if_else(agegp == '5-17y', 'gp3', 'gp4')))) %>%
  dplyr::group_by(yearc, age, stg) %>%
  dplyr::summarise(n = sum(n)) %>%
  dplyr::ungroup()

#age-specific IPD cases
ipdAfit7 <-
  ipd7 %>%
  tidyr::pivot_wider(., id_cols = yearc, names_from = c(stg, age), values_from = n) %>%
  dplyr::ungroup()  %>%
  dplyr::select(yearc, V_gp1, V_gp2, V_gp3, V_gp4, F_gp1, F_gp2, F_gp3, F_gp4, N_gp1, N_gp2, N_gp3, N_gp4) %>%
  tidyr::replace_na(list(V_gp1=0, V_gp2=0, V_gp3=0, V_gp4=0, F_gp1=0, F_gp2=0, F_gp3=0, F_gp4=0, N_gp1=0, N_gp2=0, N_gp3=0, N_gp4=0))

ipdAfit7A <-
  ipd7A %>%
  tidyr::pivot_wider(., id_cols = yearc, names_from = c(stg, age), values_from = n) %>%
  dplyr::ungroup()  %>%
  dplyr::select(yearc, V_gp1, V_gp2, V_gp3, V_gp4, F_gp1, F_gp2, F_gp3, F_gp4, N_gp1, N_gp2, N_gp3, N_gp4) %>%
  tidyr::replace_na(list(V_gp1=0, V_gp2=0, V_gp3=0, V_gp4=0, F_gp1=0, F_gp2=0, F_gp3=0, F_gp4=0, N_gp1=0, N_gp2=0, N_gp3=0, N_gp4=0))

ipdAfit13 <-
  ipd13 %>%
  tidyr::pivot_wider(., id_cols = yearc, names_from = c(stg, age), values_from = n) %>%
  dplyr::ungroup()  %>%
  dplyr::select(yearc, V_gp1, V_gp2, V_gp3, V_gp4, F_gp1, F_gp2, F_gp3, F_gp4, N_gp1, N_gp2, N_gp3, N_gp4) %>%
  tidyr::replace_na(list(V_gp1=0, V_gp2=0, V_gp3=0, V_gp4=0, F_gp1=0, F_gp2=0, F_gp3=0, F_gp4=0, N_gp1=0, N_gp2=0, N_gp3=0, N_gp4=0))

ipdAfit_all <- dplyr::bind_rows(ipdAfit7, ipdAfit13)
ipdAfit_all_A <- dplyr::bind_rows(ipdAfit7, ipdAfit7A)

#total ipd cases
ipdTfit7 <-
  ipdAfit7 %>%
  dplyr::mutate(V = rowSums(.[2:5]),
                F = rowSums(.[6:9]),
                N = rowSums(.[10:13])) %>%
  dplyr::select(yearc, V, F, N)

ipdTfit13 <-
  ipdAfit13 %>%
  dplyr::mutate(V = rowSums(.[2:5]),
                F = rowSums(.[6:9]),
                N = rowSums(.[10:13])) %>%
  dplyr::select(yearc, V, F, N)

ipdTfit13 <- dplyr::bind_rows(ipdTfit7, ipdTfit13)

#ggplot pcv7 cases
ipd7 %>%
  ggplot(aes(x = yearc, y = log(n+1))) +
  geom_line() +
  facet_grid(age ~ stg) +
  theme_bw() +
  geom_vline(aes(xintercept = 2000), linetype = "dotted", linewidth = 0.5) +
  geom_vline(aes(xintercept = 2010), linetype = "dotted", linewidth = 0.5) +
  #scale_x_discrete(limits = c(1998, 2000, 2002, 2004, 2006, 2008)) +
  labs(title = "", x = "year", y = "log_ipd") +
  theme(legend.position = "none")

#ggplot pcv13 cases
ipd13 %>%
  ggplot(aes(x = yearc, y = log(n+1))) +
  geom_line() +
  facet_grid(age ~ stg) +
  theme_bw() +
  geom_vline(aes(xintercept = 2000), linetype = "dotted", linewidth = 0.5) +
  geom_vline(aes(xintercept = 2010), linetype = "dotted", linewidth = 0.5) +
  scale_x_discrete(limits = c(2010, 2012, 2014, 2016, 2018, 2020)) +
  labs(title = "", x = "year", y = "log_ipd") +
  theme(legend.position = "none")

#ggplot pcv7+pcv13
dplyr::bind_rows(ipd7,ipd13) %>%
  ggplot(aes(x = yearc, y = log(n+1))) +
  geom_line() +
  facet_grid(age ~ stg) +
  theme_bw() +
  geom_vline(aes(xintercept = 2000), linetype = "dotted", linewidth = 0.5) +
  geom_vline(aes(xintercept = 2010), linetype = "dotted", linewidth = 0.5) +
  scale_x_discrete(limits = c(1999, 2005, 2010, 2015, 2018, 2020)) +
  labs(title = "", x = "year", y = "log_ipd") +
  theme(legend.position = "none")

#====================================================================

#case to carrier ratios
ccr <-
  dplyr::bind_cols(
    ipdAfit_all %>% dplyr::filter(yearc == 1999) %>% tidyr::pivot_longer(V_gp1:N_gp4),
    data_frame(name2 = c('V_gp1', 'V_gp2', 'V_gp3', 'V_gp4', 'F_gp1', 'F_gp2', 'F_gp3', 'F_gp4', 'N_gp1', 'N_gp2', 'N_gp3', 'N_gp4'),
               value2 = c(508850.1, 2561032.1, 8229984.2, 33763076.4, 508850.1, 2561032.1, 8229984.2, 33763076.4, 508850.1, 2561032.1, 8229984.2, 33763076.4))) %>%
  dplyr::mutate(invas = value/value2)


ccry <-
  dplyr::bind_cols(
    ipdAfit_all %>% dplyr::filter(yearc == 2010) %>% tidyr::pivot_longer(V_gp1:N_gp4),
    data_frame(name2 = c('V_gp1', 'V_gp2', 'V_gp3', 'V_gp4', 'F_gp1', 'F_gp2', 'F_gp3', 'F_gp4', 'N_gp1', 'N_gp2', 'N_gp3', 'N_gp4'),
               value2 = c(13504.76, 155396.40, 1121108.66, 5679847.39, 33673.98, 386408.14, 2778238.42, 14182369.78, 626719.66, 1808404.75, 11207576.27, 55871212.99))) %>%
  dplyr::mutate(invas = value/value2)

#====================================================================

# #daily average number of physical contacts between any two specific individuals in age i and age j (Flasche et al.)
contact_US <-
  as.matrix(rbind(c(0.40, 0.06, 0.02, 0.04),
                  c(0.80, 1.88, 0.33, 0.43),
                  c(0.53, 0.74, 2.24, 0.52),
                  c(2.37, 1.87, 1.38, 1.42))
  )
