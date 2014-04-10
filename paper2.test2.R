for (i in 1:100) {
  parms[[2]] = alphaA - i * stepwise
  lvout <- ode(y = Nstar, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
  Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
  Nstar[Nstar < 10^-8] = 0  # species with biomass less than the threshold is considered to be extinct
  
  comm = jacobian.full(y = Nstar, func = lv2, parms = parms)  # Jacobian matrix in equilibrium
  lev = max(Re(eigen(comm)$values))
  
  survived = sum(Nstar > 0)  # Survived species at equillibrim state
  extinct = sum(Nstar == 0)  # Extinct species at equillibrim state
  Nstars[,i+1] = Nstar
  print(lev)
}

##### evaluate 'Disentangling nestedness from models of ecological complex'

matplot(t, lvout[, 2:10], type = "l", lwd = 1.5, col = 1)

############ Nestedness definition
source('nestedness.R')  # import the functions defining different nestednesss

A = Safariland.C.Pre  # incident matrix
A[A > 0] = 1  # transform to binary matrix
B = as.one.mode(A)  # transform to adjacency matrix
B2 = B %*% B # get the overlap matrix, which indicate the common neighbors between species
sum(B2)  # the sum of all the common neighbors between species
degrees = rowSums(B)  # the degrees (number of neighbors) of species
degrees^2  # the number of neighbors of neighbors of species
print(paste(sum(B2), sum(degrees^2)))  # 


