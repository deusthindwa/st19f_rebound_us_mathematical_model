# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

# #====================================================================
# 
# #load observed carriage data
# source(here::here("script", "01_preprocessData.R"))
# obs_traj <- df_EW
# 
# #extract parameter sets
# parset <- rio::import(here::here("results", "lhs_parsetX.csv"))
# base1 <- parset$V1
# base2 <- parset$V2
# base3 <- parset$V3
# base4 <- parset$V4
# 
# base5 <- parset$V5
# base6 <- parset$V6
# base7 <- parset$V7
# base8 <- parset$V8
# 
# base9 <- parset$V9
# base10 <- parset$V10
# base11 <- parset$V11
# base12 <- parset$V12
# 
# # compV1 <- parset$V5
# # compV2 <- parset$V6
# # compV3 <- parset$V7
# # compV4 <- parset$V8
# # compF1 <- parset$V9
# # compF2 <- parset$V10
# # compF3 <- parset$V11
# # compF4 <- parset$V12
# # compN1 <- parset$V13
# # compN2 <- parset$V14
# # compN3 <- parset$V15
# # compN4 <- parset$V16
# 
# #load model-simulations
# #pred_traj <- Filter(function(x) !is.null(x) && length(x) > 0, pred_traj)
# pred_traj <- readRDS(here::here("results", "saved_carr_dynamicsX.rds"))
# 
# #multinomial likelihood function
# calc_carr_logLik = function(obs_count_carr, model_prev_carr){
#   LL = sum(obs_count_carr * log(model_prev_carr))
#   if(any(model_prev_carr <= 0)) LL = -1000001
#   return(LL)
# }
# 
# #multinomial likelihood for carriage given V, F, N and S states
# GOF_prev <- c()
# N_ages <- 4
# 
# for (i in 1:N_sample_size) {
#   sim_traj <- pred_traj[[i]]
#   LL_prop = 0
#   for (age in 1:N_ages){
#     LL_prop = LL_prop +
#       calc_carr_logLik(as.numeric(c(obs_traj$VT.prev[age], obs_traj$F.prev[age], obs_traj$NVT.prev[age], obs_traj$Uncol[age])) * obs_traj$N.total[age],
#                        as.numeric(c(sim_traj$V.prev.model[age], sim_traj$F.prev.model[age], sim_traj$N.prev.model[age], sim_traj$S.prev.model[age])) 
#       )
#   }
#   GOF_prev[i] <- LL_prop
# }
# 
# #select 1000 posterior samples
# XX <-
#   data_frame(data.frame(GOF_prev) %>% mutate(resampled_rows = as.integer(rownames(.)))) %>%
#   dplyr::arrange(desc(GOF_prev)) %>%
#   dplyr::mutate(rownum = as.integer(rownames(.))) %>%
#   dplyr::filter(rownum <= 1000)
# 
# resampled_rows <- XX$resampled_rows
# 
# pbase1 <- base1[resampled_rows]
# pbase2 <- base2[resampled_rows]
# pbase3 <- base3[resampled_rows]
# pbase4 <- base4[resampled_rows]
# 
# pbase5 <- base5[resampled_rows]
# pbase6 <- base6[resampled_rows]
# pbase7 <- base7[resampled_rows]
# pbase8 <- base8[resampled_rows]
# 
# pbase9 <- base9[resampled_rows]
# pbase10 <- base10[resampled_rows]
# pbase11 <- base11[resampled_rows]
# pbase12 <- base12[resampled_rows]
# 
# # pcompV1 <- compV1[resampled_rows]
# # pcompV2 <- compV2[resampled_rows]
# # pcompV3 <- compV3[resampled_rows]
# # pcompV4 <- compV4[resampled_rows]
# # pcompF1 <- compF1[resampled_rows]
# # pcompF2 <- compF2[resampled_rows]
# # pcompF3 <- compF3[resampled_rows]
# # pcompF4 <- compF4[resampled_rows]
# # pcompN1 <- compN1[resampled_rows]
# # pcompN2 <- compN2[resampled_rows]
# # pcompN3 <- compN3[resampled_rows]
# # pcompN4 <- compN4[resampled_rows]
# 
# #create dataset for each posterior parameter estimate
# posterior <-
#   bind_rows(
#     data.frame(value = (pbase1), par = "Vρ(<1y)"),
#     data.frame(value = (pbase2), par = "Vρ(1-4y)"),
#     data.frame(value = (pbase3), par = "Vρ(5-17y)"),
#     data.frame(value = (pbase4), par = "Vρ(18+y)"),
#     
#     data.frame(value = (pbase5), par = "Fρ(<1y)"),
#     data.frame(value = (pbase6), par = "Fρ(1-4y)"),
#     data.frame(value = (pbase7), par = "Fρ(5-17y)"),
#     data.frame(value = (pbase8), par = "Fρ(18+y)"),
#     
#     data.frame(value = (pbase9), par = "Nρ(<1y)"),
#     data.frame(value = (pbase10), par = "Nρ(1-4y)"),
#     data.frame(value = (pbase11), par = "Nρ(5-17y)"),
#     data.frame(value = (pbase12), par = "Nρ(18+y)")
#     
#     # data.frame(value = (pcompV1), par = "εV"),
#     # data.frame(value = (pcompV2), par = "εV2"),
#     # data.frame(value = (pcompV3), par = "εV3"),
#     # data.frame(value = (pcompV4), par = "εV4"),
#     # data.frame(value = (pcompF1), par = "εF"),
#     # data.frame(value = (pcompF2), par = "εF2"),
#     # data.frame(value = (pcompF3), par = "εF3"),
#     # data.frame(value = (pcompF4), par = "εF4"),
#     # data.frame(value = (pcompN1), par = "εN"),
#     # data.frame(value = (pcompN2), par = "εN2"),
#     # data.frame(value = (pcompN3), par = "εN3"),
#     # data.frame(value = (pcompN4), par = "εN4")
#     )
# 
# #summary of estimated posterior parameters
# posteriorEst <- 
#   posterior %>%
#   dplyr::group_by(par) %>%
#   dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
#                    l = round(quantile(value, 0.025), digits = 3),
#                    u = round(quantile(value, 0.975), digits = 3), na.rm = T)
# 
# #plot posterior parameter distributions
# # "εN",#"εN2","εN3","εN4",
# # "εF",#"εF2","εF3","εF4", 
# # "εV",#"εV2","εV3","εV4", 
# posteriorPlot <-
#   posteriorEst %>%
#   dplyr::mutate(par = factor(par, levels = c("Vρ(18+y)", "Vρ(5-17y)", "Vρ(1-4y)", "Vρ(<1y)",
#                                              "Fρ(18+y)", "Fρ(5-17y)", "Fρ(1-4y)", "Fρ(<1y)",
#                                              "Nρ(18+y)", "Nρ(5-17y)", "Nρ(1-4y)", "Nρ(<1y)")),
#                 na.rm = 'Model calibration') %>%
#   dplyr::filter(!is.na(par)) %>%
#   ggplot() +
#   geom_point(aes(x = log(m), y=par), color = 'black', size = 2, stroke=1.5, shape =1) +
#   geom_errorbar(aes(x=log(m), y=par, xmin = log(l), xmax = log(u)), color = 'black', width = 0, size = 0.8, stat = "identity") +
#   guides(color = guide_legend(title="Parameter")) +
#   facet_grid(.~na.rm) +
#   theme_bw(base_size = 14, base_family = "American Typewriter") +
#   scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4.5, 4)) +
#   #geom_text(aes(x = 4.5, y = 7, label = paste0(m[1],' (',l[7],'-',u[7],')')), hjust = 0, vjust = 0.5, size = 4, family = "serif", color = "gray10") +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 10)) +
#   labs(title = "(b)", x = "Log_estimate", y = "") +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
#   theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
#   theme(legend.position = 'none')
# 
# rio::export(posteriorEst, file = here("output", "prevalence_estimates.csv"))
# 
# #find the simulations corresponding to re-samples
# post_traj <- list()
# for (i in 1:6){
#   j =  resampled_rows[i]
#   post_traj[[i]] <- pred_traj[[j]]
# }
# 
# #create dataset for fitted prevalence
# post_traj <- 
#   dplyr::bind_rows(
#     dplyr::bind_rows(post_traj, .id = "run") %>%
#       dplyr::group_by(agegp) %>%
#       dplyr::summarise(Est = round(quantile(V.prev.model, 0.500), digits = 3),
#                        Low = round(quantile(V.prev.model, 0.025), digits = 3),
#                        High = round(quantile(V.prev.model, 0.975), digits = 3)) %>%
#       dplyr::mutate(stg = "V"),
#     
#     dplyr::bind_rows(post_traj, .id = "run") %>%
#       dplyr::group_by(agegp) %>%
#       dplyr::summarise(Est = round(quantile(F.prev.model, 0.500), digits = 3),
#                        Low = round(quantile(F.prev.model, 0.025), digits = 3),
#                        High = round(quantile(F.prev.model, 0.975), digits = 3)) %>%
#       dplyr::mutate(stg = "F"),
#     
#     dplyr::bind_rows(post_traj, .id = "run") %>%
#       dplyr::group_by(agegp) %>%
#       dplyr::summarise(Est = round(quantile(N.prev.model, 0.500), digits = 3),
#                        Low = round(quantile(N.prev.model, 0.025), digits = 3),
#                        High = round(quantile(N.prev.model, 0.975), digits = 3)) %>%
#       dplyr::mutate(stg = "N")) %>%
#   dplyr::mutate(carrtype = "Fitted")
#   
# #combine fitted and observed prevalence for ggplot
# post_trajPlot <-
#   dplyr::bind_rows(post_traj, df_EWo) %>%
#   dplyr::mutate(stg = factor(stg, levels = c('V', 'F', 'N')),
#                 agegp = factor(agegp, levels = c('<1y', '1-4y', '5-17y', '18+y'))) %>%
#   ggplot() +
#   geom_point(aes(x = agegp, y=Est, colour = carrtype), size = 2, stroke=1.5, shape =1, position = position_dodge(width = 0.2)) +
#   geom_errorbar(aes(x=agegp, y=Est,  ymin = Low, ymax = High, color=carrtype) , width = 0, size = 0.8, position = position_dodge(width = 0.2), stat = "identity") +
#   scale_color_manual(values = c('#F8766D', 'gray40')) +
#   guides(color = guide_legend(title="")) +
#   labs(title = "(c)", x = "Age group", y = "Carriage prevalence") +
#   facet_grid(. ~ stg) +
#   theme_bw(base_size = 14, base_family = "American Typewriter") +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
#   theme(legend.position = c(0.85,0.85))
# 
# post_trajPlot
# 
# ggsave(here::here("output", "prevalence_obsfit.png"),
#        plot = post_trajPlot,
#        width = 12, height = 6, unit = "in", dpi = 300)




#====================================================================

#load observed carriage data
source(here::here("script", "carriage", "01_preprocessData.R"))
obs_traj <- df_EW

#extract parameters
parset <- rio::import(here::here("results", "carriage", "lhs_parset.csv"))
base1 <- parset$V1
base2 <- parset$V2
base3 <- parset$V3
base4 <- parset$V4
compV <- parset$V5
compF <- parset$V6
compN <- parset$V7

#load model simulated data
pred_traj <- readRDS(here::here("results", "carriage", "saved_carr_dynamics.rds"))
#pred_traj <- Filter(function(x) !is.null(x) && length(x) > 0, pred_traj)

#multinomial likelihood function
calc_carr_logLik = function(obs_count_carr, model_prev_carr){
  LL = sum(obs_count_carr * log(model_prev_carr)) 
  if(any(model_prev_carr <= 0)) LL = -1000001
  return(LL)
}

#multinomial likelihood for carriage given V, F, N and S states
GOF_prev <- c()
N_sample_size <- 1000
N_ages <- 4

for (i in 1:N_sample_size) {
  sim_traj <- pred_traj[[i]]
  LL_prop=0
  for (age in 1:N_ages){
    LL_prop = LL_prop +
      calc_carr_logLik(as.numeric(c(obs_traj$VT.prev[age], obs_traj$F.prev[age], obs_traj$NVT.prev[age], obs_traj$Uncol[age])) * obs_traj$N.total[age],
                       as.numeric(c(sim_traj$V.prev.model[age], sim_traj$F.prev.model[age], sim_traj$N.prev.model[age], sim_traj$S.prev.model[age])) 
      )
  }
  GOF_prev[i] <- LL_prop
}

#select 1000 posterior samples
XX <- 
  data_frame(data.frame(GOF_prev) %>% mutate(resampled_rows = as.integer(rownames(.)))) %>%
  dplyr::arrange(desc(GOF_prev)) %>%
  dplyr::mutate(rownum = as.integer(rownames(.))) %>%
  dplyr::filter(rownum <= 5000)

resampled_rows <- XX$resampled_rows

#normalized weights to sum to 1 
#GOF_prev <- GOF_prev[!is.na(GOF_prev)]
# log_weights <- GOF_prev - max(GOF_prev)
# weights <- exp(log_weights)
# weights <- weights / sum(weights)
# n_resamples <- N_sample_size/100 # keep 1-10 ratio

#re-sample from the proposal distribution using the calculated weights
#resampled_rows <- sample.int(seq_along(compV), size = n_resamples, replace = TRUE, prob = weights)
pbase1 <- base1[resampled_rows]
pbase2 <- base2[resampled_rows]
pbase3 <- base3[resampled_rows]
pbase4 <- base4[resampled_rows]
pcompV <- compV[resampled_rows]
pcompF <- compF[resampled_rows]
pcompN <- compN[resampled_rows]

#create dataset for each posterior parameter estimate
posterior <-
  bind_rows(
    data.frame(value = (pbase1), par = "ρ(<1y)"),
    data.frame(value = (pbase2), par = "ρ(1-4y)"),
    data.frame(value = (pbase3), par = "ρ(5-17y)"),
    data.frame(value = (pbase4), par = "ρ(18+y)"),
    data.frame(value = (pcompV), par = "εV"),
    data.frame(value = (pcompF), par = "εF"),
    data.frame(value = (pcompN), par = "εN"))

#summary of estimated posterior parameters
posteriorEst <- 
  posterior %>%
  dplyr::group_by(par) %>%
  dplyr::summarise(m = round(quantile(value, 0.500), digits = 3),
                   l = round(quantile(value, 0.025), digits = 3),
                   u = round(quantile(value, 0.975), digits = 3), na.rm=T)

#plot posterior parameter distributions
posteriorPlot <-
  posteriorEst %>%
  dplyr::mutate(par = factor(par, levels = c("εN", "εF", "εV", "ρ(18+y)", "ρ(5-17y)", "ρ(1-4y)", "ρ(<1y)")),
                na.rm = 'Model calibration') %>%
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

rio::export(posteriorEst, file = here("output", "carriage", "posteriorEst.csv"))

#find the simulations corresponding to re-samples
#resampled_rows <- c(130527, 31147, 88050, 103213, 179302, 174341)
post_traj <- list()
for (i in 1:6){
  j =  resampled_rows[i]
  post_traj[[i]] <- pred_traj[[j]]
}

#create dataset for fitted prevalence
post_traj <- 
  dplyr::bind_rows(
    dplyr::bind_rows(post_traj, .id = "run") %>%
      dplyr::group_by(agegp) %>%
      dplyr::summarise(Est = round(quantile(V.prev.model, 0.500), digits = 3),
                       Low = round(quantile(V.prev.model, 0.025), digits = 3),
                       High = round(quantile(V.prev.model, 0.975), digits = 3)) %>%
      dplyr::mutate(stg = "V"),
    
    dplyr::bind_rows(post_traj, .id = "run") %>%
      dplyr::group_by(agegp) %>%
      dplyr::summarise(Est = round(quantile(F.prev.model, 0.500), digits = 3),
                       Low = round(quantile(F.prev.model, 0.025), digits = 3),
                       High = round(quantile(F.prev.model, 0.975), digits = 3)) %>%
      dplyr::mutate(stg = "F"),
    
    dplyr::bind_rows(post_traj, .id = "run") %>%
      dplyr::group_by(agegp) %>%
      dplyr::summarise(Est = round(quantile(N.prev.model, 0.500), digits = 3),
                       Low = round(quantile(N.prev.model, 0.025), digits = 3),
                       High = round(quantile(N.prev.model, 0.975), digits = 3)) %>%
      dplyr::mutate(stg = "N")) %>%
  dplyr::mutate(carrtype = "Fitted")

#combine fitted and observed prevalence for ggplot
post_trajPlot <-
  dplyr::bind_rows(post_traj, df_EWo) %>%
  dplyr::mutate(stg = factor(stg, levels = c('V', 'F', 'N')),
                agegp = factor(agegp, levels = c('<1y', '1-4y', '5-17y', '18+y'))) %>%
  ggplot() +
  geom_point(aes(x = agegp, y=Est, colour = carrtype), size = 2, stroke=1.5, shape =1, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(x=agegp, y=Est,  ymin = Low, ymax = High, color=carrtype) , width = 0, size = 0.8, position = position_dodge(width = 0.2), stat = "identity") +
  scale_color_manual(values = c('#F8766D', 'gray40')) +
  guides(color = guide_legend(title="")) +
  labs(title = "(c)", x = "Age group", y = "Carriage prevalence") +
  facet_grid(. ~ stg) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = c(0.85,0.85))

((modelstru/posteriorPlot) | post_trajPlot | plot_layout(ncol = 2, width = c(1,2)))

ggsave(here::here("output", "carriage", "obs_fit_prevalence.png"),
       plot = ((modelstru/posteriorPlot) | post_trajPlot | plot_layout(ncol = 2, width = c(1,2))),
       width = 18, height = 10, unit = "in", dpi = 300)
