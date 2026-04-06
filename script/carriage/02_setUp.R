# by Deus Thindwa
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#total population in each  age group
pop_EW <- as.numeric(c(629200, 3009200, 10131900, 39274600))
names(pop_EW) <- c('agegp1', 'agegp2', 'agegp3', 'agegp4')

#leave blank if constant population where births == deaths
birth_EW <- matrix(NA, nrow=1000, ncol=4)

#daily average number of contacts between individuals in age i and j
contact_EW <- contact_EW

#plot daily average contact matrix
#heatmap(contact_EW/sum(diag(contact_EW)), Rowv = NA, Colv = NA, scale = 'none')

#could replace this with vector of actual age names
N_ages <- length(pop_EW)
agenames <- paste0('agegp', 1:N_ages)

#define states of transmission compartmental model
StateNames <- c('S', 'V', 'F', 'N', 'NV', 'NF', 'VF')
N_states <- length(StateNames)

#N age groups x K states
yinit.matrix <- array(NA, dim = c(N_ages, length(StateNames)))

#assign row names and column names
dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

#initializes population with with 1 infected person in each infection state (V,N,F,M)
yinit.matrix[,'S'] = c(pop_EW[1:4])
yinit.matrix['agegp1','V'] = yinit.matrix['agegp1','F'] = yinit.matrix['agegp1','N'] = 1
yinit.matrix[, c('NV')] = 0
yinit.matrix[, c('NF')] = 0
yinit.matrix[, c('VF')] = 0
yinit.matrix[,'S'] = c(pop_EW[1:4]-3)
yinit.matrix[,'V'] = yinit.matrix[,'F'] = yinit.matrix[,'N'] = 1

#vectorize the states for ODE input. Vectorize the ynit matrix
yinit.vector <- as.vector(yinit.matrix)

#create array that has the labels by age and state, and use this to name the yinit.vector
name.array <- array(NA, dim = dim(yinit.matrix))
for(i in 1:dim(name.array)[1]){ # to number of rows of yinit.matrix
  for(j in 1:dim(name.array)[2]){ # to number of columns for yinit.matrix
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i], # loop through the row names from agegrp1
                             dimnames(yinit.matrix)[[2]][j]) # loop through the column names from state M
  }
}
name.vector <- as.vector(name.array) # vectorise the name.array
names(yinit.vector) <- name.vector # assign labels to yinit.vector
