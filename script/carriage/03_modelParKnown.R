# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#-------------------------------------------------------------------------------

#daily duration of carriage (r) assumed similar by ST-group but different by age-group
durV = c(71.42857, 34.48276, 17.85714, 16.94915)
durF = c(71.42857, 34.48276, 17.85714, 16.94915)
durN = c(71.42857, 34.48276, 17.85714, 16.94915)

#probability of transmission per contact (per capita) (beta)
#this is assumed for whole carriage duration hence needs divided by carriage duration
#currently assumed similar by ST-group and by age-group but will be varied and estimated
#baseV = c(53.1, 2.66, 0.54, 0.05)
#baseF = c(53.1, 2.66, 0.54, 0.05)
#baseN = c(53.1, 2.66, 0.54, 0.05)

#daily average number of physical contacts between any two specific individuals in age i and age j
#frequency dependent transmission is assumed
contact_EW = contact_EW

#susceptibility of resident ST-group to be infected by another ST-group relative to uncolonized individuals
#currently assumed similar by ST-group and by age-group but will be varied and estimated
compV = 0.5
compF = 0.5
compN = 0.5

#annual duration of being in age class
agewidth= c(1, 3, 13, 83)

#states and age classes
yinit.matrix = yinit.matrix

#number of age groups
N_ages = N_ages

#save parameters in a list
parms <- list(
  durV = durV,
  durF = durF,
  durN = durN,

  # baseV = baseV,
  # baseF = baseF,
  # baseN = baseN,
  
  contact_EW = contact_EW,

  compV = compV,
  compF = compF,
  compN = compN,
  
  agewidth = agewidth,
  yinit.matrix = yinit.matrix,
  N_ages=N_ages,
  time.step = 'year'
)
