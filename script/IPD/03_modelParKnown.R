# by Deus Thindwa and Dan Weinberger
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#daily duration of carriage (r) assumed similar by ST-group but different by age-group (should the clearance vary by serogroup as)
durV = c(71.42857, 34.48276, 17.85714, 16.94915)
durF = c(71.42857, 34.48276, 17.85714, 16.94915)
durN = c(71.42857, 34.48276, 17.85714, 16.94915)

#probability of transmission per contact (per capita) (beta)
#this is assumed for whole carriage duration hence needs divided by carriage duration
#assumed similar by ST-group but different by age-group as estimated from carriage data in England
baseV = c(42.512, 1.5, 0.25, 0.25)
baseF = c(42.512, 1.5, 0.25, 0.25)
baseN = c(42.512, 1.5, 0.25, 0.25)

# In density-dependent transmission, the contact rate (c) depends on the population density. 
# In frequent-dependent transmission, the contact rate (c) does not depend on the population density. (we assumed frequency-dependent)
contact_US = contact_US

#susceptibility of resident ST-group to be co-colonised by another ST-group relative to uncolonized individuals
#currently assumed similar by ST-group and by age-group as estimated in English carriage data
# compV = 0.5
# compF = 0.5
# compN = 0.5

#vaccination parameters 
#delta1 is %VE against V serotypes AND delta2 is %VE against F serotypes
#omega1 is VE waning against V serotypes AND omega2 is VE waning against F serotypes
# delta1 = 0.60
# delta2 = 0.60
# omega1 = 8
# omega2 = 8

#number of age groups
N_ages = N_ages

#case carrier ratios, proportion of carriage to IPD per 100,000 population (observed IPD cases/predicted carriage incidence)
invV = c(ccr$invas[1], ccr$invas[2], ccr$invas[3], ccr$invas[4])
invF = c(ccr$invas[5], ccr$invas[6], ccr$invas[7], ccr$invas[8])
invN = c(ccr$invas[9], ccr$invas[10], ccr$invas[11], ccr$invas[12])

#annual duration of being in age class (assume life expectancy is 77y in US)
agewidth = c(1, 4, 13, 59) #59

#states and age classes
yinit.matrix = yinit.matrix

#US age population
pop_US = pop_US

#burn-in (this is added to indicate burn-in period)
t0 = 11

#save parameters in a list
parms <- list(
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
  # 
  # delta1 = delta1,
  # delta2 = delta2,
  # omega1 = 8,
  # omega2 = 8,
  
  invV = invV,
  invF = invF,
  invN = invN,

  agewidth = agewidth,
  yinit.matrix = yinit.matrix,
  pop_US = pop_US,
  N_ages = N_ages,
  time.step = 'year',
  t0 = t0
  )
