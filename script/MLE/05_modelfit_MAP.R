# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#clear R environment
rm(list = ls())

#load R packages
source(here::here("script", "MLE", "00_main.R"))

#load data to fit model to
source(here::here("script", "MLE", "01_preprocessData.R"))
obs_ipd <- ipdAfit7 %>% slice_head(., n=11) %>% dplyr::select(everything(), -yearc)

#set up model including births, contacts, population
source(here::here("script", "MLE", "02_setUp.R"))

#set up fixed parameter values
source(here::here("script", "MLE", "03_modelParKnown.R"))

#run model function
source(here::here("script", "MLE", "04_transModel.R"))


#set all fixed parameters
my_parmset <- list(
  durV = durV,
  durF = durF,
  durN = durN,
  
  baseV = baseV,
  baseF = baseF,
  baseN = baseN,
  
  contact_US = contact_US,
  
  # compV = compV,
  # compF = compF,
  # compN = compN,
  
  invV = invV,
  invF = invF,
  invN = invN,
  
  pop_US = pop_US,
  obs_pop_1999_2009 = obs_pop_1999_2009,
  obs_pop_2010_2019 = obs_pop_2010_2019,
  obs_pop_1999_2019 = obs_pop_1999_2019,
  
  agewidth = agewidth,
  yinit.matrix = yinit.matrix,
  N_ages = N_ages,
  time.step = 'year',
  t0 = t0
)

#time step is in years
start_time = 1
tmax =  500
my_times <- seq(start_time, tmax, by = 1)



fitmodel <- function(parameters, dat) {

#parameters to estimate
  
  #vaccination parameters
  delta1 <- parameters[1]
  delta2 <- parameters[2]
  omega1 <- parameters[3]
  omega2 <- parameters[4]
  compV <-  parameters[5]
  compF <-  parameters[6]
  compN <-  parameters[7]
  replacex <- parameters[(8:11)]

  #we need this (double arrow assignment) for the function and model to recognize these parameters
  delta1 <<- delta1
  delta2 <<- delta2
  baseN <<- baseN
  omega1 <<- omega1
  omega2 <<- omega2
  compV <<- compV
  compF <<- compF
  compN <<- compN
  replacex <<- replacex

  #run transmission model with initial conditions & time steps defined above, and parameter values from function call
  results <- ode(y = yinit.vector, 
                 t = my_times,
                 func = simple_model,
                 parms = c(my_parmset,
                           list(delta1 = delta1,
                                delta2 = delta2,
                                omega1 = omega1,
                                omega2 = omega2,
                                compV = compV,
                                compF = compF,
                                compN = compN,
                                replacex = replacex
                                )
                           )
                 )

  #extract time varying force of infection
  lambdaDS <-
    data.frame(results) %>%
    dplyr::select(X57:X68)
  
  #extract states results
  results <-
    data.frame(results) %>%
    dplyr::select(everything(), -X57:-X68)
  
  #evaluation period excludes burn-in period
  t0 <- my_parmset$t0
  results <- tail(results, t0)
  lambdaDS <- tail(lambdaDS, t0)
  
  steadyStates <- results[,-1] #remove time variable
  
  #unvaccinated states
  S <-  steadyStates[,grep('\\b(S)\\b',  colnames(steadyStates))]
  V <-  steadyStates[,grep('\\b(V)\\b',  colnames(steadyStates))]
  F <-  steadyStates[,grep('\\b(F)\\b',  colnames(steadyStates))]
  N <-  steadyStates[,grep('\\b(N)\\b',  colnames(steadyStates))]
  NV <- steadyStates[,grep('\\b(NV)\\b', colnames(steadyStates))]
  NF <- steadyStates[,grep('\\b(NF)\\b', colnames(steadyStates))]
  VF <- steadyStates[,grep('\\b(VF)\\b', colnames(steadyStates))]
  
  #PCV7 vaccinated states
  Sx <-  steadyStates[,grep('\\b(Sx)\\b',  colnames(steadyStates))]
  Vx <-  steadyStates[,grep('\\b(Vx)\\b',  colnames(steadyStates))]
  Fx <-  steadyStates[,grep('\\b(Fx)\\b',  colnames(steadyStates))]
  Nx <-  steadyStates[,grep('\\b(Nx)\\b',  colnames(steadyStates))]
  NVx <- steadyStates[,grep('\\b(NVx)\\b', colnames(steadyStates))]
  NFx <- steadyStates[,grep('\\b(NFx)\\b', colnames(steadyStates))]
  VFx <- steadyStates[,grep('\\b(VFx)\\b', colnames(steadyStates))]
  
  if (my_parmset$time.step == 'year') {
    length.step = 365.25} else if (my_parmset$time.step == 'month') {
      length.step = 30.44} else {length.step = 7} #weekly
  
  #number of age groups
  N_ages = my_parmset$N_ages
  
  #clearance rates
  clearV = 1/(my_parmset$durV/length.step)
  clearF = 1/(my_parmset$durF/length.step)
  clearN = 1/(my_parmset$durN/length.step)
  
  #force of infection at steady states
  lambdaV = lambdaF = lambdaN = matrix(0, nrow = t0, ncol = my_parmset$N_ages)
  lambdaV <- lambdaDS[,1:4]
  lambdaF <- lambdaDS[,5:8]
  lambdaN <- lambdaDS[,9:12]
  
  # #vaccination parameters
  # delta1 = parms$delta1 #VE against PCV7 serotypes
  # delta2 = parms$delta2 #VE against 19F
  
  VEd = 0.95 #VE against PCV7 disease or 19F disease
  VEv = 1-(1-VEd)/(1-delta1) #VE against progression to PCV7 disease
  VEf = 1-(1-VEd)/(1-delta2) #VE against progression to 19F disease
  
  #no of new carriers in all carrier states
  ipdV = ipdF = ipdN = matrix(0, nrow = t0, ncol = N_ages)
  for (i in 1:N_ages){
    for (t in 1) {
      ipdV[t,i] = (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i])
      ipdF[t,i] = (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i])
      ipdN[t,i] = (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i])
    }
  }
  
  #initial infections across age groups
  ccr$value2 <- c(ipdV[1,], ipdF[1,], ipdN[1,])
  ccr <- ccr %>% dplyr::mutate(invas = value/value2)
  source(here::here("script", "MLE", "03_modelParKnown.R"))
  
  #case:carrier ratio (= observed incidence/steady state incidence before vaccination)
  #this is calculated in '03_modelParKnown' file
  invV = my_parmset$invV
  invF = my_parmset$invF
  invN = my_parmset$invN
  
  #adjusted for case-carrier ratios before vaccination
  for (i in 1:N_ages){
    for (t in 1) {
      ipdV[t,i] = invV[i] * ipdV[t,i]
      ipdF[t,i] = invF[i] * ipdF[t,i]
      ipdN[t,i] = invN[i] * ipdN[t,i]
    }
  }
  
  #adjusted for case-carrier ratios after vaccination
  for (i in 1:N_ages){
    for (t in 2:t0) {
      ipdV[t,i] = invV[i] * (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i]) * (1-VEv))
      ipdF[t,i] = invF[i] * (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i]) * (1-VEf))
      ipdN[t,i] = invN[i] * (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i]) ) #* replacex[i]
    }
  }
  
  #total cases by serotype group and age groups
  ipdStateAge <- cbind(ipdV, ipdF, ipdN)
  colnames(ipdStateAge) <- c('V_gp1','V_gp2','V_gp3','V_gp4',
                             'F_gp1','F_gp2','F_gp3','F_gp4',
                             'N_gp1','N_gp2','N_gp3','N_gp4')
  
  #convert entire numeric dataset to integer for easy poisson distribution
  sim_traj <<- as_tibble(lapply(as_tibble(ipdStateAge), as.integer))
  
  #read onserved data to fit the model to
  obs_traj <<- dat
  obs_pop <<- my_parmset$obs_pop_1999_2009
  

  #likelihood under poisson distribution
  LL_pois <- 0
  
  for (k in 1:t0){
    LL_pois = LL_pois +
      sum(dpois(x = obs_traj$V_gp1[k], lambda = sim_traj$V_gp1[k], log = TRUE),
          dpois(x = obs_traj$V_gp2[k], lambda = sim_traj$V_gp2[k], log = TRUE),
          dpois(x = obs_traj$V_gp3[k], lambda = sim_traj$V_gp3[k], log = TRUE),
          dpois(x = obs_traj$V_gp4[k], lambda = sim_traj$V_gp4[k], log = TRUE),
          
          dpois(x = obs_traj$F_gp1[k], lambda = sim_traj$F_gp1[k], log = TRUE),
          dpois(x = obs_traj$F_gp2[k], lambda = sim_traj$F_gp2[k], log = TRUE),
          dpois(x = obs_traj$F_gp3[k], lambda = sim_traj$F_gp3[k], log = TRUE),
          dpois(x = obs_traj$F_gp4[k], lambda = sim_traj$F_gp4[k], log = TRUE),
          
          dpois(x = obs_traj$N_gp1[k], lambda = sim_traj$N_gp1[k], log = TRUE),
          dpois(x = obs_traj$N_gp2[k], lambda = sim_traj$N_gp2[k], log = TRUE),
          dpois(x = obs_traj$N_gp3[k], lambda = sim_traj$N_gp3[k], log = TRUE),
          dpois(x = obs_traj$N_gp4[k], lambda = sim_traj$N_gp4[k], log = TRUE)
      )
  }
  
  # for (k in 1:t0){
  #   LL_pois = LL_pois +
  #     sum(dpois(x = obs_traj$V_gp1[k], lambda = sim_traj$V_gp1[k]*obs_pop$V_pop1[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$V_gp2[k], lambda = sim_traj$V_gp2[k]*obs_pop$V_pop2[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$V_gp3[k], lambda = sim_traj$V_gp3[k]*obs_pop$V_pop3[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$V_gp4[k], lambda = sim_traj$V_gp4[k]*obs_pop$V_pop4[k]/1e6, log = TRUE),
  #         
  #         dpois(x = obs_traj$F_gp1[k], lambda = sim_traj$F_gp1[k]*obs_pop$F_pop1[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$F_gp2[k], lambda = sim_traj$F_gp2[k]*obs_pop$F_pop2[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$F_gp3[k], lambda = sim_traj$F_gp3[k]*obs_pop$F_pop3[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$F_gp4[k], lambda = sim_traj$F_gp4[k]*obs_pop$F_pop4[k]/1e6, log = TRUE),
  #         
  #         dpois(x = obs_traj$N_gp1[k], lambda = sim_traj$N_gp1[k]*obs_pop$N_pop1[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$N_gp2[k], lambda = sim_traj$N_gp2[k]*obs_pop$N_pop2[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$N_gp3[k], lambda = sim_traj$N_gp3[k]*obs_pop$N_pop3[k]/1e6, log = TRUE),
  #         dpois(x = obs_traj$N_gp4[k], lambda = sim_traj$N_gp4[k]*obs_pop$N_pop4[k]/1e6, log = TRUE)
  #     )
  # }
  
  #prior distributions
  delta1_prior <- dunif(x = delta1, min = 0.6, max = 0.9, log = TRUE) #PCV VE against V
  delta2_prior <- dunif(x = delta2, min = 0.6, max = 0.9, log = TRUE) #PCV VE against F
  omega1_prior <- dunif(x = omega1, min = 5, max = 20, log = TRUE) #PCV waning rate against V
  omega2_prior <- dunif(x = omega2, min = 5, max = 20, log = TRUE) #PCV waning rate against F
  compV_prior  <- dunif(x = compV, min = 0, max = 1, log = TRUE) #V relative risk to compete F, N
  compF_prior  <- dunif(x = compF, min = 0, max = 1, log = TRUE) #F relative risk to compete V, N
  compN_prior  <- dunif(x = compN, min = 0, max = 1, log = TRUE) #N relative risk to compete V, F
  replacex_prior1 <- dunif(x = replacex[1], min = 1, max = 4, log = TRUE) #N relative risk to compete V, F
  replacex_prior2 <- dunif(x = replacex[2], min = 1, max = 10, log = TRUE) #N relative risk to compete V, F
  replacex_prior3 <- dunif(x = replacex[3], min = 1, max = 10, log = TRUE) #N relative risk to compete V, F
  replacex_prior4 <- dunif(x = replacex[4], min = 100, max = 150, log = TRUE) #N relative risk to compete V, F
  
  #total log-likelihood (sum up because of log) = likelihood and priors
  LL_total <- LL_pois + delta1_prior + delta2_prior + omega1_prior + omega2_prior + compV_prior + compF_prior + compN_prior + replacex_prior1 + replacex_prior2 + replacex_prior3 + replacex_prior4
  
  return(LL_total)
  
}


#you can minimize the negative log likelihood in order to actually perform the maximum likelihood estimate of the function you are testing.
#use Optim function for maximum a posteriori (calibrate outputs of the transmission model to observed status)
#starting values for parameters
fitLL <- optim(par = c(0.7, 0.7, 12, 12, 0.5, 0.5, 0.5, 1.5, 5, 5, 120),
               fn = fitmodel,
               dat = obs_ipd, #("dat" argument is passed to the function specified in fn)
               control = list(fnscale = -1, maxit = 30000000, trace = TRUE, REPORT = 500))


# fitLL <- minqa::bobyqa(
#   par = c(0.7, 0.7, 12, 12, 0.5, 0.5, 0.5, 5, 5, 5, 120), 
#   fn = fitmodel,
#   dat = obs_ipd,
#   control = list(maxfun=100000)
# )

#check if the convergence = 0. 
#if not run for several times to make sure the function did not stuck at local minimum
fitLL
fitLL$par
#exp(fitLL$par)

# var_cov_matrix <- solve(fitLL$hessian, tol = .Machine$double.eps.)
# std_errors <- sqrt(diag(var_cov_matrix))
# fitLL$par - 1.96 * std_errors
# fitLL$par + 1.96 * std_errors

#====================================================================
#====================================================================


#extract fitted parameters
parms <- c(my_parmset,
           list(delta1 = (fitLL$par[1]),
           delta2 = (fitLL$par[2]),
           omega1 = (fitLL$par[3]),
           omega2 = (fitLL$par[4]),
           compV  = (fitLL$par[5]),
           compF  = (fitLL$par[6]),
           compN  = (fitLL$par[7]),
           replacex = c(fitLL$par[8], fitLL$par[9], fitLL$par[10], fitLL$par[11]))
           )

#check calibration results
results <- ode(y = yinit.vector, 
               t = my_times,
               func = simple_model,
               parms = parms)

#extract time varying force of infection
lambdaDS <-
  data.frame(results) %>%
  dplyr::select(X57:X68)

#extract states results
results <-
  data.frame(results) %>%
  dplyr::select(everything(), -X57:-X68)

#evaluation period excludes burn-in period
t0 <- parms$t0
results <- tail(results, t0)
lambdaDS <- tail(lambdaDS, t0)

steadyStates <- results[,-1] #remove time variable

#unvaccinated states
S <-  steadyStates[,grep('\\b(S)\\b',  colnames(steadyStates))]
V <-  steadyStates[,grep('\\b(V)\\b',  colnames(steadyStates))]
F <-  steadyStates[,grep('\\b(F)\\b',  colnames(steadyStates))]
N <-  steadyStates[,grep('\\b(N)\\b',  colnames(steadyStates))]
NV <- steadyStates[,grep('\\b(NV)\\b', colnames(steadyStates))]
NF <- steadyStates[,grep('\\b(NF)\\b', colnames(steadyStates))]
VF <- steadyStates[,grep('\\b(VF)\\b', colnames(steadyStates))]

#PCV7 vaccinated states
Sx <-  steadyStates[,grep('\\b(Sx)\\b',  colnames(steadyStates))]
Vx <-  steadyStates[,grep('\\b(Vx)\\b',  colnames(steadyStates))]
Fx <-  steadyStates[,grep('\\b(Fx)\\b',  colnames(steadyStates))]
Nx <-  steadyStates[,grep('\\b(Nx)\\b',  colnames(steadyStates))]
NVx <- steadyStates[,grep('\\b(NVx)\\b', colnames(steadyStates))]
NFx <- steadyStates[,grep('\\b(NFx)\\b', colnames(steadyStates))]
VFx <- steadyStates[,grep('\\b(VFx)\\b', colnames(steadyStates))]

if (parms$time.step == 'year') {
  length.step = 365.25} else if (parms$time.step == 'month') {
    length.step = 30.44} else {length.step = 7} #weekly

#number of age groups
N_ages = parms$N_ages

#clearance rates
clearV = 1/(parms$durV/length.step)
clearF = 1/(parms$durF/length.step)
clearN = 1/(parms$durN/length.step)

#force of infection at steady states
lambdaV = lambdaF = lambdaN = matrix(0, nrow = t0, ncol = parms$N_ages)
lambdaV <- lambdaDS[,1:4]
lambdaF <- lambdaDS[,5:8]
lambdaN <- lambdaDS[,9:12]

#vaccination parameters
delta1 = parms$delta1 #VE against PCV7 serotypes
delta2 = parms$delta2 #VE against 19F

VEd = 0.95 #VE against PCV7 disease or 19F disease
VEv = 1-(1-VEd)/(1-delta1) #VE against progression to PCV7 disease
VEf = 1-(1-VEd)/(1-delta2) #VE against progression to 19F disease

#no of new carriers in all carrier states
ipdV = ipdF = ipdN = matrix(0, nrow = t0, ncol = N_ages)
for (i in 1:N_ages){
  for (t in 1) {
    ipdV[t,i] = (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i])
    ipdF[t,i] = (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i])
    ipdN[t,i] = (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i])
  }
}

#initial infections across age groups
ccr$value2 <- c(ipdV[1,], ipdF[1,], ipdN[1,])
ccr <- ccr %>% dplyr::mutate(invas = value/value2)
source(here::here("script", "MLE", "03_modelParKnown.R"))

#case:carrier ratio (= observed incidence/steady state incidence before vaccination)
#this is calculated in '03_modelParKnown' file
invV = parms$invV
invF = parms$invF
invN = parms$invN

#adjusted for case-carrier ratios before vaccination
for (i in 1:N_ages){
  for (t in 1) {
    ipdV[t,i] = invV[i] * ipdV[t,i]
    ipdF[t,i] = invF[i] * ipdF[t,i]
    ipdN[t,i] = invN[i] * ipdN[t,i]
  }
}

#adjusted for case-carrier ratios after vaccination
for (i in 1:N_ages){
  for (t in 2:t0) {
    ipdV[t,i] = invV[i] * (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i]) * (1-VEv))
    ipdF[t,i] = invF[i] * (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i]) * (1-VEf))
    ipdN[t,i] = invN[i] * (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i]) * replacex[i])
  }
}

#total cases by serotype group and age groups
ipdStateAge <- cbind(ipdV, ipdF, ipdN)
colnames(ipdStateAge) <- c('V_gp1','V_gp2','V_gp3','V_gp4',
                           'F_gp1','F_gp2','F_gp3','F_gp4',
                           'N_gp1','N_gp2','N_gp3','N_gp4')

#convert entire numeric dataset to integer for easy poisson distribution
posteriorDS <- 
  as_tibble(lapply(as_tibble(ipdStateAge), as.integer)) %>% dplyr::mutate(yearc = (seq.int(1999, 2009, by=1))) %>%
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
                cat = 'Predicted') %>%
  dplyr::ungroup() %>%
  dplyr::rename('m'='value')


observedDS_xt <-
ipdAfit7 %>%
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

