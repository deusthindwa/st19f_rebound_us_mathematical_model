# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#load observed carriage data
source(here::here("script", "IPD", "01_preprocessData.R"))
obs_traj <- ipdAfit7 %>% slice_head(., n=12)

#extract parameters
parset <- rio::import(here::here("results", "IPD", "lhs_parset.csv"))
delta1 <- parset$V1
delta2 <- parset$V2
omega1 <- parset$V3
omega2 <- parset$V4

rep1 <- parset$V8
rep2 <- parset$V9
rep3 <- parset$V10
rep4 <- parset$V11

#load model simulated data
#pred_traj <- Filter(function(x) !is.null(x) && length(x) > 0, pred_traj)
pred_traj <- readRDS(here::here("results", "IPD", "saved_ipd_dynamics.rds"))

#likelihood under poisson distribution
GOF_cases <- c()
N_sample_size <- 500

for (i in 1:N_sample_size) {
  sim_traj <- pred_traj[[i]]
  LL_pois = 0
  for (j in 1:t0){
    LL_pois = LL_pois +
      sum(dpois(x = obs_traj$V_gp1[j], lambda = sim_traj$V_gp1[j], log = TRUE),
          dpois(x = obs_traj$V_gp2[j], lambda = sim_traj$V_gp2[j], log = TRUE),
          dpois(x = obs_traj$V_gp3[j], lambda = sim_traj$V_gp3[j], log = TRUE),
          dpois(x = obs_traj$V_gp4[j], lambda = sim_traj$V_gp4[j], log = TRUE),
          
          dpois(x = obs_traj$F_gp1[j], lambda = sim_traj$F_gp1[j], log = TRUE),
          dpois(x = obs_traj$F_gp2[j], lambda = sim_traj$F_gp2[j], log = TRUE),
          dpois(x = obs_traj$F_gp3[j], lambda = sim_traj$F_gp3[j], log = TRUE),
          dpois(x = obs_traj$F_gp4[j], lambda = sim_traj$F_gp4[j], log = TRUE),

          dpois(x = obs_traj$N_gp1[j], lambda = sim_traj$N_gp1[j], log = TRUE),
          dpois(x = obs_traj$N_gp2[j], lambda = sim_traj$N_gp2[j], log = TRUE),
          dpois(x = obs_traj$N_gp3[j], lambda = sim_traj$N_gp3[j], log = TRUE),
          dpois(x = obs_traj$N_gp4[j], lambda = sim_traj$N_gp4[j], log = TRUE)
          )
  }
  GOF_cases[i] <- LL_pois
}

#select 1000 posterior samples with minimum log-Likelihood
n_resamples <- N_sample_size/5
XX <-
  data_frame(data.frame(GOF_cases) %>% mutate(resampled_rows = as.integer(rownames(.)))) %>%
  dplyr::arrange(desc(GOF_cases)) %>%
  dplyr::mutate(rownum = as.integer(rownames(.))) %>%
  dplyr::filter(rownum <= n_resamples)

resampled_rows <- XX$resampled_rows

#assign samples to parameters
pdelta1 <- delta1[resampled_rows]
pdelta2 <- delta2[resampled_rows]
pomega1 <- omega1[resampled_rows]
pomega2 <- omega2[resampled_rows]

prep1 <- rep1[resampled_rows]
prep2 <- rep2[resampled_rows]
prep3 <- rep3[resampled_rows]
prep4 <- rep4[resampled_rows]

#create dataset for each parameter
posterior <-
  bind_rows(
    data.frame(value = (pdelta1), par = 'δ1'),
    data.frame(value = (pdelta2), par = 'δ2'),
    data.frame(value = (pomega1), par = 'ω1'),
    data.frame(value = (pomega2), par = 'ω2'),
    
    data.frame(value = (prep1), par = 'r1'),
    data.frame(value = (prep2), par = 'r2'),
    data.frame(value = (prep3), par = 'r3'),
    data.frame(value = (prep4), par = 'r4')
    )

#summary of estimated parameters
posteriorEst <- 
  posterior %>%
  dplyr::group_by(par) %>%
  dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
                   l = round(quantile(value, 0.025), digits = 3),
                   u = round(quantile(value, 0.975), digits = 3), na.rm=T)

posteriorPlot <-
  posteriorEst %>%
  dplyr::filter(!par %in% c('r1, r2, r3, r4')) %>%
  dplyr::mutate(par = factor(par, levels = c('δ1', 'δ2', 'ω1', 'ω2')),
                na.rm = 'Parameter estimation') %>%
  ggplot() +
  geom_point(aes(x = log(m), y=par), color = 'black', size = 2, stroke=1.5, shape =1) +
  geom_errorbar(aes(x=log(m), y=par, xmin = log(l), xmax = log(u)), color = 'black', width = 0, size = 0.8, stat = "identity") +
  guides(color = guide_legend(title="Parameter")) +
  facet_grid(.~na.rm) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10)) +
  labs(title = "(b)", x = "Log_estimate", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'none')

rio::export(posteriorEst, file = here("output", "IPD", "baseModel_estimates.csv"))

#find the simulations corresponding to resamples
posterior_traj <- list()
for (i in 1:n_resamples){
  j = resampled_rows[i]
  posterior_traj[[i]] <- pred_traj[[j]] %>% dplyr::mutate(yearc = c(1999:2009))
}

posteriorDS <-
  dplyr::bind_rows(posterior_traj, .id = "run") %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
                   l = round(quantile(value, 0.025), digits = 3),
                   u = round(quantile(value, 0.975), digits = 3), na.rm=T) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(yearc = as.integer(yearc),
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y',
                                if_else(agegp == 'gp2', '1-4y',
                                        if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                yearc = as.integer(yearc),
                cat = 'Fitted')

observedDS <-
  obs_traj %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::mutate(value = as.integer(value)) %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::mutate(l = stats::poisson.test(value, conf.level = 0.95 )$conf.int[1],
                u = stats::poisson.test(value, conf.level = 0.95 )$conf.int[2],
                na.rm = TRUE,
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y',
                                       if_else(agegp == 'gp2', '1-4y',
                                               if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                cat = 'Observed') %>%
  dplyr::ungroup() %>%
  dplyr::rename('m'='value')

observedDS_xt <-
  ipdAfit %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::mutate(value = as.integer(value)) %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::mutate(l = stats::poisson.test(value, conf.level = 0.95 )$conf.int[1],
                u = stats::poisson.test(value, conf.level = 0.95 )$conf.int[2],
                na.rm = TRUE,
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y',
                                       if_else(agegp == 'gp2', '1-4y',
                                               if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                cat = 'Observed') %>%
  dplyr::ungroup() %>%
  dplyr::rename('m'='value')

post_trajPlot <-
  ggplot() +
  geom_point(data = observedDS_xt, aes(x = yearc, y = m), shape = 1, stroke = 1.5) + 
  geom_line(data = posteriorDS, aes(x = yearc, y = m), color = '#00BFC4') +
  scale_x_continuous(breaks = seq(1999, 2019, 4), limits = c(1998, 2019)) +
  geom_ribbon(data = posteriorDS, aes(x = yearc, y = m, ymin = l, ymax = u), fill = '#00BFC4', alpha = 0.2, stat = "identity") +
  geom_errorbar(data = observedDS_xt, aes(x = yearc, y = m, ymin = l, ymax = u), color = 'black', width = 0, size = 0.6, ) +
  geom_vline(aes(x = yearc, y = m), xintercept = 2009.5, linetype = 'dashed') +
  facet_wrap(stg ~ agegp, scales = 'free') + 
  theme_bw(base_size = 14, base_family = "American Typewriter") + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10)) +
  labs(title = "", x = "Year", y = "Reported number of IPD cases") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'none')

post_trajPlot

ggsave(here::here("output", "IPD", "baseModel_obsfit.png"),
       plot = post_trajPlot,
       width = 18, height = 10, unit = "in", dpi = 300)
