# by Deus & Dan & Ginny
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

#import observed ipd data
source(here::here("script", "IPD", "scenarios", "01_preprocessData_scenario.R"))
obs_traj <- ipdAfit_all %>% slice_head(., n = 21)
obs_trajA <- ipdAfit_all_A %>% slice_head(., n = 21)

for (vcase in 1:10) {
  
scenario = vcase

#import LHS parameter set data
if (scenario == 1){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_projection.csv"))
  
} else if (scenario == 2){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_efficacy.csv"))
  
} else if (scenario == 3){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_duration.csv"))
  
} else if (scenario == 4){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_competition.csv"))
  
} else if (scenario == 5){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_modefit_efficacy.csv"))
  
} else if (scenario == 6){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_modefit_duration.csv"))
  
} else if (scenario == 7){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_modefit_competition.csv"))
  
} else if (scenario == 8){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_efficacy_reduced.csv"))
  
} else if (scenario == 9){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_efficacy_constant.csv"))
  
} else if (scenario == 10){
  parset <- rio::import(here::here("results", "IPD", "lhs_parset_project_pcv7serotypes.csv"))
  
} else{
  print("No parameter set available")
}

#extra individual priors
prior_delta1  = parset$V1
prior_delta2  = parset$V2
prior_delta3  = parset$V3
prior_delta4  = parset$V4
prior_omega1  = parset$V5
prior_omega2  = parset$V6
prior_omega3  = parset$V7
prior_omega4  = parset$V8
prior_compVx  = parset$V9
prior_compFx  = parset$V10
prior_compNx  = parset$V11
prior_compVy  = parset$V12
prior_compFy  = parset$V13
prior_compNy  = parset$V14

#import model predictions data
if (scenario == 1){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_projection.rds"))
  
} else if (scenario == 2){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_efficacy.rds"))
  
} else if (scenario == 3){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_duration..rds"))
  
} else if (scenario == 4){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_competition.rds"))
  
} else if (scenario == 5){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_modefit_efficacy.rds"))
  
} else if (scenario == 6){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_modefit_duration.rds"))
  
} else if (scenario == 7){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_modefit_competition.rds"))
  
} else if (scenario == 8){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_efficacy_reduced.rds"))
  
} else if (scenario == 9){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_efficacy_constant.rds"))
  
} else if (scenario == 10){
  pred_traj <- base::readRDS(here::here("results", "IPD", "saved_ipdRun_project_pcv7serotypes.rds"))
  
} else{
  print("No ipdRun to import")
}

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

#select posterior samples with minimum log-Likelihood
n_resamples <- N_sample_size/5
XX <-
  data_frame(data.frame(GOF_cases) %>% mutate(resampled_rows = as.integer(rownames(.)))) %>%
  dplyr::arrange(desc(GOF_cases)) %>%
  dplyr::mutate(rownum = as.integer(rownames(.))) %>%
  dplyr::filter(rownum <= n_resamples)

#logLikelihood for scenarios of model fit to the data
if (scenario == 5) { #model fit of reduced efficacy
  rio::export(data_frame(data.frame(XX$GOF_cases)), here::here("results", "IPD", "logLik_efficacy.csv"))
  
} else if (scenario == 6){
  rio::export(data_frame(data.frame(XX$GOF_cases)), here::here("results", "IPD", "logLik_duration.csv"))
  
} else if (scenario == 7){
  rio::export(data_frame(data.frame(XX$GOF_cases)), here::here("results", "IPD", "logLik_competition.csv"))
  
} else{
  print('No logLik dataset')
}

#sampled rows
resampled_rows <- XX$resampled_rows

#assign samples to parameters
prior_delta1  = prior_delta1[resampled_rows]
prior_delta2  = prior_delta2[resampled_rows]
prior_delta3  = prior_delta3[resampled_rows]
prior_delta4  = prior_delta4[resampled_rows]
prior_omega1  = prior_omega1[resampled_rows]
prior_omega2  = prior_omega2[resampled_rows]
prior_omega3  = prior_omega3[resampled_rows]
prior_omega4  = prior_omega4[resampled_rows]
prior_compVx  = prior_compVx[resampled_rows]
prior_compFx  = prior_compFx[resampled_rows]
prior_compNx  = prior_compNx[resampled_rows]
prior_compVy  = prior_compVy[resampled_rows]
prior_compFy  = prior_compFy[resampled_rows]
prior_compNy  = prior_compNy[resampled_rows]

#create dataset for each parameter
posterior <-
  bind_rows(
    data.frame(value = (prior_delta1), par = 'δ1'),
    data.frame(value = (prior_delta2), par = 'δ2'),
    data.frame(value = (prior_delta3), par = 'δ3'),
    data.frame(value = (prior_delta4), par = 'δ4'),
    data.frame(value = (prior_omega1), par = 'w1'),
    data.frame(value = (prior_omega2), par = 'w2'),
    data.frame(value = (prior_omega3), par = 'w3'),
    data.frame(value = (prior_omega4), par = 'w4'),
    data.frame(value = (prior_compVx), par = 'eVx'),
    data.frame(value = (prior_compFx), par = 'eFx'),
    data.frame(value = (prior_compNx), par = 'eNx'),
    data.frame(value = (prior_compVy), par = 'eVy'),
    data.frame(value = (prior_compFy), par = 'eFy'),
    data.frame(value = (prior_compNy), par = 'eNy')
  )

#summary of estimated parameters
posteriorEst <- 
  posterior %>%
  dplyr::group_by(par) %>%
  dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
                   l = round(quantile(value, 0.025), digits = 3),
                   u = round(quantile(value, 0.975), digits = 3), 
                   na.rm = TRUE)

if (scenario == 1){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_projection.csv"))
  
} else if (scenario == 2){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_project_efficacy.csv"))
  
} else if (scenario == 3){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_project_duration.csv"))
  
} else if (scenario == 4){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_project_competition.csv"))
  
} else if (scenario == 5){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_efficacy.csv"))
  
} else if (scenario == 6){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_duration.csv"))
  
} else if (scenario == 7){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_competition.csv"))
  
} else if (scenario == 8){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_efficacy_reduced.csv"))
  
} else if (scenario == 9){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_efficacy_constant.csv"))
  
} else if (scenario == 10){
  rio::export(posteriorEst, file = here("output", "IPD", "estimates_modefit_pcv7serotypes.csv"))
  
} else{
  print("No ipdRun to import")
}

#find the simulations corresponding to resamples
posterior_traj <- list()
for (i in 1:n_resamples){
  j = resampled_rows[i]
  posterior_traj[[i]] <- pred_traj[[j]] %>% dplyr::mutate(yearc = c(1999:2019))
}

posteriorDS <-
  dplyr::bind_rows(posterior_traj, .id = "run") %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
                   l = round(quantile(value, 0.025), digits = 3),
                   u = round(quantile(value, 0.975), digits = 3), 
                   na.rm = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(yearc = as.integer(yearc),
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y', if_else(agegp == 'gp2', '1-4y', if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                yearc = as.integer(yearc),
                cat = 'Fitted') %>%
  dplyr::mutate(timex = if_else(yearc <=2009, 'prepcv13', 'postpcv13'))

observedDS <-
  obs_traj %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::mutate(value = as.integer(value)) %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::mutate(l = stats::poisson.test(value, conf.level = 0.95 )$conf.int[1],
                u = stats::poisson.test(value, conf.level = 0.95 )$conf.int[2],
                na.rm = TRUE,
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y', if_else(agegp == 'gp2', '1-4y', if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                cat = 'Observed') %>%
  dplyr::ungroup() %>%
  dplyr::rename('m'='value')

post_trajPlot <-
  ggplot() +
  geom_point(data = observedDS, aes(x = yearc, y = m), shape = 1, stroke =1.5) + 
  geom_line(data = posteriorDS, aes(x = yearc, y = m, color = timex, group = timex), size = 0.8) +
  geom_ribbon(data = posteriorDS, aes(x = yearc, y = m, ymin = l, ymax = u, group = timex, fill = timex), alpha = 0.2, stat = "identity") +
  geom_errorbar(data = observedDS, aes(x = yearc, y = m, ymin = l, ymax = u), color = 'black', width = 0, size = 0.6, ) +
  scale_x_continuous(breaks  = seq(1999, 2019, by = 5),  limits = c(1999, 2019),  
                     labels  = function(x){gsub(pattern = "1999", replacement = "1998/9", x = x)}) +
  geom_vline(aes(x = yearc, y = m), xintercept = 2009.5, linetype = 'dashed') +
  facet_wrap(stg ~ agegp, scales = 'free') + 
  theme_bw(base_size = 14, base_family = "American Typewriter") + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10)) +
  labs(title = "", x = "Year", y = "Reported number of IPD cases") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'none')

#---

observedDS <-
  obs_trajA %>%
  tidyr::pivot_longer(cols = V_gp1:N_gp4, names_to = "stg") %>%
  dplyr::mutate(value = as.integer(value)) %>%
  dplyr::group_by(yearc, stg) %>%
  dplyr::mutate(l = stats::poisson.test(value, conf.level = 0.95 )$conf.int[1],
                u = stats::poisson.test(value, conf.level = 0.95 )$conf.int[2],
                na.rm = TRUE,
                agegp = stringr::str_sub(stg, 3, nchar(stg)),
                agegp = factor(if_else(agegp == 'gp1', '<1y', if_else(agegp == 'gp2', '1-4y', if_else(agegp == 'gp3', '5-17y', '18+y'))), levels = c('<1y','1-4y','5-17y','18+y')),
                stg = factor(stringr::str_sub(stg, 1, 1), levels = c('V','F','N')),
                cat = 'Observed') %>%
  dplyr::ungroup() %>%
  dplyr::rename('m'='value')

post_trajPlotA <-
  ggplot() +
  geom_point(data = observedDS, aes(x = yearc, y = m), shape = 1, stroke =1.5) + 
  geom_line(data = posteriorDS, aes(x = yearc, y = m, group = timex, color = timex), size = 0.8) +
  geom_ribbon(data = posteriorDS, aes(x = yearc, y = m, ymin = l, ymax = u, group = timex, fill = timex), alpha = 0.15, stat = "identity") +
  geom_errorbar(data = observedDS, aes(x = yearc, y = m, ymin = l, ymax = u), color = 'black', width = 0, size = 0.6, ) +
  scale_x_continuous(breaks  = seq(1999, 2019, by = 5),  limits = c(1999, 2019),  
                     labels  = function(x){gsub(pattern = "1999", replacement = "1998/9", x = x)}) +
  geom_vline(aes(x = yearc, y = m), xintercept = 2009.5, linetype = 'dashed') +
  facet_wrap(stg ~ agegp, scales = 'free') + 
  theme_bw(base_size = 14, base_family = "American Typewriter") + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10)) +
  labs(title = "", x = "Year", y = "Reported number of IPD cases") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'none')

#---

#view plot
post_trajPlot

#save all plots
if (scenario == 1){
  ggsave(here::here("output", "IPD", "projections.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 2){
  ggsave(here::here("output", "IPD", "project_efficacy.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 3){
  ggsave(here::here("output", "IPD", "project_duration.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 4){
  ggsave(here::here("output", "IPD", "project_competition.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 5){
  ggsave(here::here("output", "IPD", "modefit_efficacy.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 6){
  ggsave(here::here("output", "IPD", "modefit_duration.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 7){
  ggsave(here::here("output", "IPD", "modefit_competition.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 8){
  ggsave(here::here("output", "IPD", "project_efficacy_reduced.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 9){
  ggsave(here::here("output", "IPD", "project_efficacy_constant.png"), plot = post_trajPlot, width = 20, height = 10, unit = "in", dpi = 300)
  
} else if (scenario == 10){
  ggsave(here::here("output", "IPD", "project_efficacy_pcv7serotypes.png"), plot = post_trajPlotA, width = 20, height = 10, unit = "in", dpi = 300)
  
}

else{
  print("No ipdRun to import")
}

}


#====================================================================

Model_logLik <-
  dplyr::bind_cols(
    logLik_efficacy <- rio::import(here::here('results', "IPD", 'logLik_efficacy.csv')) %>% dplyr::rename('ll_efficacy' = 'XX.GOF_cases'),
    logLik_duration <- rio::import(here::here('results', "IPD", 'logLik_duration.csv'))  %>% dplyr::rename('ll_duration' = 'XX.GOF_cases'),
    logLik_competition <- rio::import(here::here('results', "IPD", 'logLik_competition.csv'))  %>% dplyr::rename('ll_competition' = 'XX.GOF_cases')) %>%
  dplyr::summarise(ll_eff = mean(ll_efficacy),
                ll_dur = mean(ll_duration),
                ll_com = mean(ll_competition))

#compute AIC and BIC
# AIC= −2⋅logL + 2k
# BIC= −2⋅log(L)+ k⋅log(n)
k = 14
n = 252
Model_logLik <-
  Model_logLik %>%
  dplyr::mutate(aic_eff = -2*Model_logLik$ll_eff + 2*k,
                aic_dur = -2*Model_logLik$ll_dur + 2*k,
                aic_com = -2*Model_logLik$ll_com + 2*k,
                
                bic_eff = -2*Model_logLik$ll_eff + k*log(n),
                bic_dur = -2*Model_logLik$ll_dur + k*log(n),
                bic_com = -2*Model_logLik$ll_com + k*log(n))
