# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

fitmodel <- function(parameters, dat) {   # takes the parameter values and dataset as inputs
  
  # Run transmission model with initial conditions & time steps defined above, and parameter values from function call
  results <- ode(y = yinit.vector, 
                 t = my_times,
                 func = simple_model,
                 parms = my_parmset
  )
  
  #model extraction and prevalence
  res_burned <- results[nrow(results), , drop = FALSE]
  S <- res_burned[,grep('\\b(S)\\b', colnames(res_burned))]
  V <- res_burned[,grep('\\b(V)\\b', colnames(res_burned))]
  F <- res_burned[,grep('\\b(F)\\b', colnames(res_burned))]
  N <- res_burned[,grep('\\b(N)\\b', colnames(res_burned))]
  NV <- res_burned[,grep('\\b(NV)\\b', colnames(res_burned))]
  NF <- res_burned[,grep('\\b(NF)\\b', colnames(res_burned))]
  VF <- res_burned[,grep('\\b(VF)\\b', colnames(res_burned))]
  
  #total age-specific carriers after burnin in each obs data available infection states
  #dual carriage are equally shared between states
  S_fit <- S
  V_fit <- V + 0.5*(NV+VF)
  F_fit <- F + 0.5*(NF+VF)
  N_fit <- N + 0.5*(NV+NF)
  
  #total age-specific infections after burnin
  totInf <- (V_fit + F_fit + N_fit)
  names(totInf) <- c('agegp1 Inf', 'agegp2 Inf', 'agegp3 Inf', 'agegp4 Inf')
  
  #total age-specific pop after burnin
  totPop <- (V_fit + F_fit + N_fit + S_fit)
  names(totPop) <- c('agegp1 Tot', 'agegp2 Tot', 'agegp3 Tot', 'agegp4 Tot')
  
  #simulated prevalence
  Ps = S_fit/totPop
  Pv = V_fit/totPop
  Pf = F_fit/totPop
  Pn = N_fit/totPop
  
  df_EWm <-
    data_frame(agegp = c("<1y", "1-4y", "5-17y", "18+y"),
               V.prev.model = Pv,
               F.prev.model = Pf,
               N.prev.model = Pn,
               S.prev.model = Ps)
  
  return(df_EWm)
}

