# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission dynamics and impact of vaccination policies
# 01/09/2024

#====================================================================

fitmodel <- function(parameters, dat) {   # takes the parameter values and dataset as inputs
  
  #run transmission model with initial conditions & time steps defined above, and parameter values from function call
  results <- ode(y = yinit.vector, 
                 t = my_times,
                 func = simple_model,
                 parms = my_parmset
                 )
  
  #extract time varying force of infection
  lambdaDS <-
    data.frame(results) %>%
    dplyr::select(X85:X96)

  #extract states results
  results <-
    data.frame(results) %>%
    dplyr::select(everything(), -X85:-X96)

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
  
  #PCV13 vaccinated states
  Sy <-  steadyStates[,grep('\\b(Sy)\\b',  colnames(steadyStates))]
  Vy <-  steadyStates[,grep('\\b(Vy)\\b',  colnames(steadyStates))]
  Fy <-  steadyStates[,grep('\\b(Fy)\\b',  colnames(steadyStates))]
  Ny <-  steadyStates[,grep('\\b(Ny)\\b',  colnames(steadyStates))]
  NVy <- steadyStates[,grep('\\b(NVy)\\b', colnames(steadyStates))]
  NFy <- steadyStates[,grep('\\b(NFy)\\b', colnames(steadyStates))]
  VFy <- steadyStates[,grep('\\b(VFy)\\b', colnames(steadyStates))]
  
  if (my_parmset$time.step == 'year') {
    length.step = 365.25} else if (my_parmset$time.step == 'month') {
      length.step = 30.44} else {length.step = 7} #weekly
  
  #number of age groups
  N_ages = my_parmset$N_ages
  
  #clearance rates
  clearV = 1/(my_parmset$durV/length.step)
  clearF = 1/(my_parmset$durF/length.step)
  clearN = 1/(my_parmset$durN/length.step)
  
  lambdaV = lambdaF = lambdaN = matrix(0, nrow = t0, ncol = my_parmset$N_ages)
  lambdaV <- lambdaDS[,1:4]
  lambdaF <- lambdaDS[,5:8]
  lambdaN <- lambdaDS[,9:12]
  
  #vaccination parameters
  delta1 = my_parmset$delta1 #PCV7 VE against V serotypes
  delta2 = my_parmset$delta2 #PCV7 VE against F
  delta3 = my_parmset$delta3 #PCV13 VE against V serotypes
  delta4 = my_parmset$delta4 #PCV13 VE against F serotype
  
  VEd7 = 0.95 #VE against PCV7 disease or 19F disease
  VEd13 = 0.95 #VE against PCV13 disease or 19F disease
  VE1 = 1-(1-VEd7)/(1-delta1) #PCV7 VE against progression to V disease
  VE2 = 1-(1-VEd7)/(1-delta2) #PCV7 VE against progression to F disease
  VE3 = 1-(1-VEd13)/(1-delta3) #PCV13 VE against progression to V disease
  VE4 = 1-(1-VEd13)/(1-delta4) #PCV13 VE against progression to F disease
  replacex = my_parmset$replacex
  replacey = my_parmset$replacey
  
  #no of new carriers in all carrier states in 1998/99
  ipdV = ipdF = ipdN = matrix(0, nrow = t0, ncol = N_ages)
  for (i in 1:N_ages){
    for (t in 1) {
      ipdV[t,i] = (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i])
      ipdF[t,i] = (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i])
      ipdN[t,i] = (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i])
    }
  }
  
  #case:carrier ratio (= observed incidence/steady state incidence before vaccination)
  invV = my_parmset$invV
  invF = my_parmset$invF
  invN = my_parmset$invN
  
  invVy = my_parmset$invVy
  invFy = my_parmset$invFy
  invNy = my_parmset$invNy
  
  #adjusted for case-carrier ratios before vaccination
  for (i in 1:N_ages){
    for (t in 1) {
      ipdV[t,i] = invV[i] * ipdV[t,i]
      ipdF[t,i] = invF[i] * ipdF[t,i]
      ipdN[t,i] = invN[i] * ipdN[t,i]
    }
  }
  
  #adjusted for case-carrier ratios after vaccination 2000-2009, then 2010-2019
  for (i in 1:N_ages){
    
    for (t in 2:11) {
      ipdV[t,i] = invV[i] * (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i]) * (1-VE1))
      ipdF[t,i] = invF[i] * (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i]) * (1-VE2))
      ipdN[t,i] = invN[i] * (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i])*replacex[i])
    }
    
    for (t in 12:t0) {
      ipdV[t,i] = invVy[i] * (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sy[t,i]*lambdaV[t,i]*(1-delta3) + VFy[t,i]*clearF[i] + NVy[t,i]*clearN[i]) * (1-VE3))
      ipdF[t,i] = invF[i]  * (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sy[t,i]*lambdaF[t,i]*(1-delta4) + VFy[t,i]*clearV[i] + NFy[t,i]*clearN[i]) * (1-VE4))
      ipdN[t,i] = invNy[i] * (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sy[t,i]*lambdaN[t,i] + NFy[t,i]*clearF[i] + NVy[t,i]*clearV[i])*replacey[i])
    }
    
  }

  #total cases by serotype group and age groups
  ipdStateAge <- cbind(ipdV, ipdF, ipdN)
  colnames(ipdStateAge) <- c('V_gp1','V_gp2','V_gp3','V_gp4',
                             'F_gp1','F_gp2','F_gp3','F_gp4',
                             'N_gp1','N_gp2','N_gp3','N_gp4')
  
  #convert entire numeric dataset to integer for easy poisson likelihood
  ipdStateAge <<-
    as_tibble(lapply(as_tibble(ipdStateAge), as.integer))
  
  return(ipdStateAge)
}
