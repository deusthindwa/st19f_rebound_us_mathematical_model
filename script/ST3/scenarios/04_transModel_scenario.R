# by Deus Thindwa and Dan Weinberger
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

# Source function of the transmission dynamic model
simple_model <- function(t, y, parms, time.step = 'year'){
  
  #read in initial states and their names
  States <- array(y, dim = dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  #unify the time unit of parameter inputs
  if (parms$time.step == 'year') {
    period = 1
    length.step = 365.25 #days
  } else if (parms$time.step == 'month') {
    period = 12
    length.step = 30.44 #days
  } else { #if (parms$time.step == 'week'){
    period = 52.18
    length.step = 7 #days
  }
  
  #extract (unvaccinated) states for the model as vectors
  S  <-  States[,'S']
  V  <-  States[,'V']
  F  <-  States[,'F']
  N  <-  States[,'N']
  NV <-  States[,'NV']
  NF <-  States[,'NF']
  VF <-  States[,'VF']
  
  #extract (PCV7-vaccinated) states for the model as vectors
  Sx  <-  States[,'Sx']
  Vx  <-  States[,'Vx']
  Fx  <-  States[,'Fx']
  Nx  <-  States[,'Nx']
  NVx <-  States[,'NVx']
  NFx <-  States[,'NFx']
  VFx <-  States[,'VFx']
  
  #extract (PCV13-vaccinated) states for the model as vectors
  Sy  <-  States[,'Sy']
  Vy  <-  States[,'Vy']
  Fy  <-  States[,'Fy']
  Ny  <-  States[,'Ny']
  NVy <-  States[,'NVy']
  NFy <-  States[,'NFy']
  VFy <-  States[,'VFy']
  
  #number of age groups (randomly picked state S but would be the same for each state)
  N.ages <- length(S)
  
  #ageing rate
  if (parms$time.step == 'year') {
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  } else if (parms$time.step == 'month') {
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  } else { #if (parms$time.step == 'week'){
    ageing = 1/(parms$agewidth*period)
    ageout = ageing
    agein = c(1/sum(1/ageing), ageing[-N.ages])
  }
  
  #clearance rates of VT, 19F and NVT carriage (r) (fixed)
  clearV = 1/(parms$durV/length.step)
  clearF = 1/(parms$durF/length.step)
  clearN = 1/(parms$durN/length.step)
  
  #competition parameters (fixed and dynamic)
  if (t <= 491){
    compV = parms$compVx 
    compF = parms$compFx 
    compN = parms$compNx 
  } else if (t > 491 & t <= 500){
    compV = parms$compVx
    compF = parms$compFx
    compN = parms$compNx
  } else {
    compV = parms$compVy
    compF = parms$compFy
    compN = parms$compNy
  }

  #US population
  pop = parms$pop_US
  
  #vaccination program starts after a defined time
  #vaccination coverage (only individuals in age 0-1y are given one-time vaccine)
  if (t <= 491) {
    delta1 = 0
    delta2 = 0
    omega1 = 0
    omega2 = 0
    vaccx = c(0,0,0,0)
    delta3 = 0
    delta4 = 0
    omega3 = 0
    omega4 = 0
    pcv_switch = 0
  } else if (t > 491 & t <= 500){
    delta1 = parms$delta1
    delta2 = parms$delta2
    omega1 = 1/parms$omega1
    omega2 = 1/parms$omega2
    vaccx = c(0.95,0,0,0)
    delta3 = 0
    delta4 = 0
    omega3 = 0
    omega4 = 0
    pcv_switch = 0
  } else {
    delta1 = parms$delta1
    delta2 = parms$delta2
    omega1 = 0
    omega2 = 0
    vaccx = c(0.95,0,0,0)
    delta3 = parms$delta3
    delta4 = parms$delta4
    omega3 = 1/parms$omega3
    omega4 = 1/parms$omega4
    pcv_switch = 1
  }
  
  #transmission probability per contact
  baseV = parms$baseV
  baseF = parms$baseF
  baseN = parms$baseN
  
  #transmission probability per contact per unit time
  bV <- baseV/(parms$durV/length.step)
  bF <- baseF/(parms$durF/length.step)
  bN <- baseN/(parms$durN/length.step)
  
  #let q be transmission type (=0 when density or =1 when frequency-dependent)
  q = 1
  
  #annual mean number of contacts (from daily contacts)
  contact_US = parms$contact_US*length.step
  
  # Transmission probability per unit time in each age group 
  betaV <- (bV/100)/(sum(yinit.matrix)^(1-q))*contact_US
  betaF <- (bF/100)/(sum(yinit.matrix)^(1-q))*contact_US
  betaN <- (bN/100)/(sum(yinit.matrix)^(1-q))*contact_US
  
  #total number of carriers from each ST-group
  carryV <- (V+VF+NV+Vx+VFx+NVx+Vy+VFy+NVy)/sum(States)
  carryF <- (F+VF+NF+Fx+VFx+NFx+Fy+VFy+NFy)/sum(States)
  carryN <- (N+NF+NV+Nx+NFx+NVx+Ny+NFy+NVy)/sum(States)
  
  #sum across the columns of contact matrix (beta*I) to compute the force of infection
  lambdaV <- carryV %*% betaV
  lambdaF <- carryF %*% betaF
  lambdaN <- carryN %*% betaN
  
  #force of infection plus vectorize force of infection
  lambdaV <- as.vector(lambdaV)
  lambdaF <- as.vector(lambdaF)
  lambdaN <- as.vector(lambdaN)
  
  #create a matrix to record the changing variables
  dy <- matrix(NA, nrow = N_ages, ncol = ncol(States))
  colnames(dy) <- colnames(States)
  
  #ordinary differential equations (unvaccinated)
  
  dy[,'S'] <- clearV*V + clearN*N + clearF*F - lambdaV*S - lambdaN*S - lambdaF*S + omega1*Sx + omega3*Sy - ageout*S + agein*c(sum(States), S[-N.ages])*(1-vaccx)

  dy[,'V'] <- lambdaV*S + clearN*NV + clearF*VF - compV*lambdaN*V - compV*lambdaF*V - clearV*V + omega1*Vx + omega3*Vy - ageout*V + agein*c(0,V[-N.ages])*(1-vaccx)

  dy[,'F'] <- lambdaF*S + clearN*NF + clearV*VF - compF*lambdaN*F - compF*lambdaV*F - clearF*F + omega2*Fx + omega4*Fy - ageout*F + agein*c(0,F[-N.ages])*(1-vaccx)

  dy[,'N'] <- lambdaN*S + clearF*NF + clearV*NV - compN*lambdaV*N - compN*lambdaF*N - clearN*N + omega1*Nx + omega3*Ny - ageout*N + agein*c(0,N[-N.ages])*(1-vaccx)

  dy[,'NV'] <- compV*lambdaN*V + compN*lambdaV*N - clearV*NV - clearN*NV + omega1*NVx + omega3*NVy - ageout*NV + agein*c(0,NV[-N.ages])*(1-vaccx)

  dy[,'NF'] <- compN*lambdaF*N + compF*lambdaN*F - clearF*NF - clearN*NF + omega2*NFx + omega4*NFy - ageout*NF + agein*c(0,NF[-N.ages])*(1-vaccx)

  dy[,'VF'] <- compV*lambdaF*V + compF*lambdaV*F - clearV*VF - clearF*VF + omega2*VFx + omega4*VFy - ageout*VF + agein*c(0,VF[-N.ages])*(1-vaccx)

  
  #ordinary differential equations (PCV7-vaccinated)
  
  dy[,'Sx'] <- clearV*Vx + clearN*Nx + clearF*Fx - (1-delta1)*lambdaV*Sx - lambdaN*Sx - (1-delta2)*lambdaF*Sx - omega1*Sx - ageout*Sx + agein*c(0,Sx[-N.ages]) + agein*c(sum(States),S[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'Vx'] <- (1-delta1)*lambdaV*Sx + clearN*NVx + clearF*VFx - compV*lambdaN*Vx - compV*(1-delta2)*lambdaF*Vx - clearV*Vx - omega1*Vx - ageout*Vx + agein*c(0,Vx[-N.ages]) + agein*c(0,V[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'Fx'] <- (1-delta2)*lambdaF*Sx + clearN*NFx + clearV*VFx - compF*lambdaN*Fx - compF*(1-delta1)*lambdaV*Fx - clearF*Fx - omega2*Fx - ageout*Fx + agein*c(0,Fx[-N.ages]) + agein*c(0,F[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'Nx'] <- lambdaN*Sx + clearF*NFx + clearV*NVx - compN*(1-delta1)*lambdaV*Nx - compN*(1-delta2)*lambdaF*Nx - clearN*Nx - omega1*Nx - ageout*Nx + agein*c(0,Nx[-N.ages]) + agein*c(0,N[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'NVx'] <- compV*lambdaN*Vx + compN*(1-delta1)*lambdaV*Nx - clearV*NVx - clearN*NVx - omega1*NVx - ageout*NVx + agein*c(0,NVx[-N.ages]) + agein*c(0,NV[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'NFx'] <- compN*(1-delta2)*lambdaF*Nx + compF*lambdaN*Fx - clearF*NFx - clearN*NFx - omega2*NFx - ageout*NFx + agein*c(0,NFx[-N.ages]) + agein*c(0,NF[-N.ages])*vaccx*(1-pcv_switch)

  dy[,'VFx'] <- compV*(1-delta2)*lambdaF*Vx + compF*(1-delta1)*lambdaV*Fx - clearV*VFx - clearF*VFx - omega2*VFx - ageout*VFx + agein*c(0,VFx[-N.ages]) + agein*c(0,VF[-N.ages])*vaccx*(1-pcv_switch)

  
  #ordinary differential equations (PCV13-vaccinated)
  
  dy[,'Sy'] <- clearV*Vy + clearN*Ny + clearF*Fy - (1-delta3)*lambdaV*Sy - lambdaN*Sy - (1-delta4)*lambdaF*Sy - omega3*Sy - ageout*Sy + agein*c(0,Sy[-N.ages]) + agein*c(sum(States),S[-N.ages])*vaccx*pcv_switch
  
  dy[,'Vy'] <- (1-delta3)*lambdaV*Sy + clearN*NVy + clearF*VFy - compV*lambdaN*Vy - compV*(1-delta4)*lambdaF*Vy - clearV*Vy - omega3*Vy - ageout*Vy + agein*c(0,Vy[-N.ages]) + agein*c(0,V[-N.ages])*vaccx*pcv_switch
  
  dy[,'Fy'] <- (1-delta4)*lambdaF*Sy + clearN*NFy + clearV*VFy - compF*lambdaN*Fy - compF*(1-delta3)*lambdaV*Fy - clearF*Fy - omega4*Fy - ageout*Fy + agein*c(0,Fy[-N.ages]) + agein*c(0,F[-N.ages])*vaccx*pcv_switch
  
  dy[,'Ny'] <- lambdaN*Sy + clearF*NFy + clearV*NVy - compN*(1-delta3)*lambdaV*Ny - compN*(1-delta4)*lambdaF*Ny - clearN*Ny - omega3*Ny - ageout*Ny + agein*c(0,Ny[-N.ages]) + agein*c(0,N[-N.ages])*vaccx*pcv_switch
  
  dy[,'NVy'] <- compV*lambdaN*Vy + compN*(1-delta3)*lambdaV*Ny - clearV*NVy - clearN*NVy - omega3*NVy - ageout*NVy + agein*c(0,NVy[-N.ages]) + agein*c(0,NV[-N.ages])*vaccx*pcv_switch
  
  dy[,'NFy'] <- compN*(1-delta4)*lambdaF*Ny + compF*lambdaN*Fy - clearF*NFy - clearN*NFy - omega4*NFy - ageout*NFy + agein*c(0,NFy[-N.ages]) + agein*c(0,NF[-N.ages])*vaccx*pcv_switch
  
  dy[,'VFy'] <- compV*(1-delta4)*lambdaF*Vy + compF*(1-delta3)*lambdaV*Fy - clearV*VFy - clearF*VFy - omega4*VFy - ageout*VFy + agein*c(0,VFy[-N.ages]) + agein*c(0,VF[-N.ages])*vaccx*pcv_switch
  
  
  derivs <- as.vector(dy)
  res <- list(derivs, lambdaV, lambdaF, lambdaN)
  return(res)
}

#time step is in years
start_time = 1
tmax =  510
my_times <- seq(start_time, tmax, by = 1)

#=============================================================================

# # #EXAMINE THE TRAJECTORY WITHOUT FITTING FIRST
# 
# #using the ODE function to get simulated transmission dynamic results
# results <- ode(y = yinit.vector,
#                t = my_times,
#                func = simple_model,
#                parms = parms)
# plot(results)
# 
# #extract time varying force of infection
# lambdaDS <-
#   data.frame(results) %>%
#   dplyr::select(X85:X96)
# 
# #extract states results
# results <-
#   data.frame(results) %>%
#   dplyr::select(everything(), -X85:-X96)
# 
# #evaluation period excludes burn-in period
# my_parmset <- parms
# t0 <- my_parmset$t0
# results <- tail(results, t0)
# lambdaDS <- tail(lambdaDS, t0)
# 
# steadyStates <- results[,-1] #remove time variable
# 
# #unvaccinated states
# S <-  steadyStates[,grep('\\b(S)\\b',  colnames(steadyStates))]
# V <-  steadyStates[,grep('\\b(V)\\b',  colnames(steadyStates))]
# F <-  steadyStates[,grep('\\b(F)\\b',  colnames(steadyStates))]
# N <-  steadyStates[,grep('\\b(N)\\b',  colnames(steadyStates))]
# NV <- steadyStates[,grep('\\b(NV)\\b', colnames(steadyStates))]
# NF <- steadyStates[,grep('\\b(NF)\\b', colnames(steadyStates))]
# VF <- steadyStates[,grep('\\b(VF)\\b', colnames(steadyStates))]
# 
# #PCV7 vaccinated states
# Sx <-  steadyStates[,grep('\\b(Sx)\\b',  colnames(steadyStates))]
# Vx <-  steadyStates[,grep('\\b(Vx)\\b',  colnames(steadyStates))]
# Fx <-  steadyStates[,grep('\\b(Fx)\\b',  colnames(steadyStates))]
# Nx <-  steadyStates[,grep('\\b(Nx)\\b',  colnames(steadyStates))]
# NVx <- steadyStates[,grep('\\b(NVx)\\b', colnames(steadyStates))]
# NFx <- steadyStates[,grep('\\b(NFx)\\b', colnames(steadyStates))]
# VFx <- steadyStates[,grep('\\b(VFx)\\b', colnames(steadyStates))]
# 
# #PCV13 vaccinated states
# Sy <-  steadyStates[,grep('\\b(Sy)\\b',  colnames(steadyStates))]
# Vy <-  steadyStates[,grep('\\b(Vy)\\b',  colnames(steadyStates))]
# Fy <-  steadyStates[,grep('\\b(Fy)\\b',  colnames(steadyStates))]
# Ny <-  steadyStates[,grep('\\b(Ny)\\b',  colnames(steadyStates))]
# NVy <- steadyStates[,grep('\\b(NVy)\\b', colnames(steadyStates))]
# NFy <- steadyStates[,grep('\\b(NFy)\\b', colnames(steadyStates))]
# VFy <- steadyStates[,grep('\\b(VFy)\\b', colnames(steadyStates))]
# 
# if (my_parmset$time.step == 'year') {
#   length.step = 365.25
#   } else if (my_parmset$time.step == 'month') {
#     length.step = 30.44
#     } else {length.step = 7} #weekly
# 
# #number of age groups
# N_ages = my_parmset$N_ages
# 
# #clearance rates
# clearV = 1/(my_parmset$durV/length.step)
# clearF = 1/(my_parmset$durF/length.step)
# clearN = 1/(my_parmset$durN/length.step)
# 
# lambdaV = lambdaF = lambdaN = matrix(0, nrow = t0, ncol = my_parmset$N_ages)
# lambdaV <- lambdaDS[,1:4]
# lambdaF <- lambdaDS[,5:8]
# lambdaN <- lambdaDS[,9:12]
# 
# #vaccination parameters
# delta1 = my_parmset$delta1 #PCV7 VE against PCV7 serotypes
# delta2 = my_parmset$delta2 #PCV7 VE against 19F
# delta3 = my_parmset$delta3 #PCV13 VE against V serotypes
# delta4 = my_parmset$delta4 #PCV13 VE against F serotype
# 
# VEd7 = 0.95 #VE against PCV7 disease or 19F disease
# VEd13 = 0.95 #VE against PCV13 disease or 19F disease
# VE1 = 1-((1-VEd7)/(1-delta1)) #PCV7 VE against progression to V disease
# VE2 = 1-((1-VEd7)/(1-delta2)) #PCV7 VE against progression to F disease
# VE3 = 1-((1-VEd13)/(1-delta3)) #PCV13 VE against progression to V disease
# VE4 = 1-((1-VEd13)/(1-delta4)) #PCV13 VE against progression to F disease
# 
# #no of new carriers in all carrier states
# ipdV = ipdF = ipdN = ipdVx = ipdFx = ipdNx = matrix(0, nrow = t0, ncol = N_ages)
# for (i in 1:N_ages){
#   for (t in 1) {
#     ipdV[t,i] = (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i])
#     ipdF[t,i] = (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i])
#     ipdN[t,i] = (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i])
#   }
# }
# 
# #initial infections across age groups
# c(ipdV[1,], ipdF[1,], ipdN[1,])
# 
# #case:carrier ratio (= observed incidence/steady state incidence before vaccination)
# invV = my_parmset$invV
# invF = my_parmset$invF
# invN = my_parmset$invN
# 
# invVy = my_parmset$invVy
# invFy = my_parmset$invFy
# invNy = my_parmset$invNy
# 
# #adjusted for case-carrier ratios before vaccination
# for (i in 1:N_ages){
#   for (t in 1) {
#     ipdV[t,i] = invV[i] * ipdV[t,i]
#     ipdF[t,i] = invF[i] * ipdF[t,i]
#     ipdN[t,i] = invN[i] * ipdN[t,i]
#   }
# }
# 
# #adjusted for case-carrier ratios after vaccination
# for (i in 1:N_ages){
# 
#   for (t in 2:11) {
#     ipdV[t,i] = invV[i] * (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sx[t,i]*lambdaV[t,i]*(1-delta1) + VFx[t,i]*clearF[i] + NVx[t,i]*clearN[i]) * (1-VE1))
#     ipdF[t,i] = invF[i] * (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sx[t,i]*lambdaF[t,i]*(1-delta2) + VFx[t,i]*clearV[i] + NFx[t,i]*clearN[i]) * (1-VE2))
#     ipdN[t,i] = invN[i] * (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sx[t,i]*lambdaN[t,i] + NFx[t,i]*clearF[i] + NVx[t,i]*clearV[i]))
#   }
# 
#   for (t in 12:t0) {
#     ipdV[t,i] = (S[t,i]*lambdaV[t,i] + VF[t,i]*clearF[i] + NV[t,i]*clearN[i] + (Sy[t,i]*lambdaV[t,i]*(1-delta3) + VFy[t,i]*clearF[i] + NVy[t,i]*clearN[i]) * (1-VE3))
#     ipdF[t,i] = (S[t,i]*lambdaF[t,i] + VF[t,i]*clearV[i] + NF[t,i]*clearN[i] + (Sy[t,i]*lambdaF[t,i]*(1-delta4) + VFy[t,i]*clearV[i] + NFy[t,i]*clearN[i]) * (1-VE4))
#     ipdN[t,i] = (S[t,i]*lambdaN[t,i] + NF[t,i]*clearF[i] + NV[t,i]*clearV[i] + (Sy[t,i]*lambdaN[t,i] + NFy[t,i]*clearF[i] + NVy[t,i]*clearV[i]))
#   }
# 
# }
# 
# #initial infections across age groups
# c(ipdV[12,], ipdF[12,], ipdN[12,])
