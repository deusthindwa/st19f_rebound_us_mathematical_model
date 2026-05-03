# by Deus Thindwa and Dan Weinberger
# age-structured mathematical model for pneumococcal transmission
# 01/09/2024

#====================================================================

#total population of the US counties where ABC is in place
#pop_US <- as.numeric(c(0.02262483, 0.02649127, 0.16532220, 0.78556170)*18041585)
pop_US <- as.numeric(c(408188, 477945, 2982675, 14172778))
names(pop_US) <- c('agegp1', 'agegp2', 'agegp3', 'agegp4')

#could replace this with vector of actual age names
N_ages <- length(pop_US)
agenames <- paste0('agegp', 1:N_ages)

#define states of transmission compartmental model
StateNames <- c('S', 'V', 'F', 'N', 'NV', 'NF', 'VF', 
                'Sx', 'Vx', 'Fx', 'Nx', 'NVx', 'NFx', 'VFx',
                'Sy', 'Vy', 'Fy', 'Ny', 'NVy', 'NFy', 'VFy')
N_states <- length(StateNames)

#N age groups x K states
yinit.matrix <- array(NA, dim = c(N_ages, length(StateNames)))

# Assign row names and column names
dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

#initializes population with with 1 infected person in each infection state (V,F,N,VF,NF,NV)
yinit.matrix[, c('S')]  = c(pop_US[1:4]-6)
yinit.matrix[, c('V')]  = 1
yinit.matrix[, c('F')]  = 1
yinit.matrix[, c('N')]  = 1
yinit.matrix[, c('VF')] = 1
yinit.matrix[, c('NF')] = 1
yinit.matrix[, c('NV')] = 1

yinit.matrix[, c('Sx')]  = 0
yinit.matrix[, c('Vx')]  = 0
yinit.matrix[, c('Fx')]  = 0
yinit.matrix[, c('Nx')]  = 0
yinit.matrix[, c('VFx')] = 0
yinit.matrix[, c('NFx')] = 0
yinit.matrix[, c('NVx')] = 0

yinit.matrix[, c('Sy')]  = 0
yinit.matrix[, c('Vy')]  = 0
yinit.matrix[, c('Fy')]  = 0
yinit.matrix[, c('Ny')]  = 0
yinit.matrix[, c('VFy')] = 0
yinit.matrix[, c('NFy')] = 0
yinit.matrix[, c('NVy')] = 0

#vectorize the states for ODE input. Vectorize the ynit matrix
yinit.vector <- as.vector(yinit.matrix)

#create array that has the labels by age and state, and use this to name the yinit.vector
name.array <- array(NA, dim = dim(yinit.matrix))
for(i in 1:dim(name.array)[1]){ #to number of rows of yinit.matrix
  for(j in 1:dim(name.array)[2]){ # =to number of columns for yinit.matrix
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i], # loop through the row names from agegp1
                             dimnames(yinit.matrix)[[2]][j]) # loop through the column names from state S
  }
}
name.vector <- as.vector(name.array) # vectorise the name.array
names(yinit.vector) <- name.vector # assign labels to yinit.vector
