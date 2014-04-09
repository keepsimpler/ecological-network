library(deSolve)
library(rootSolve)

###############################################################################
#' @title Lotka-Volterra (LV) Equations of Holling type I
#' @param t time
#' @param N the variables of the LV system, a vector
#' @param parms parameters passed to LV model
#' r the intrinsic growth rate of species, a vector
#' M the interaction matrix between species
#' @return the derivation
#' @details .
#' @import deSolve
lv.1 <- function(t, N, parms) {
  r = parms[[1]]
  M = parms[[2]]
  dN <- N * (r + M %*% N)
  list(c(dN))
}

###############################################################################
#' @title Lotka-Volterra (LV) Equations of Holling type II by Bastolla et al.
#' @param t time
#' @param N the variables of the LV system, a vector
#' @param parms parameters passed to LV model
#' alphaP,alphaA the intrinsic growth rate of Plants and Animals(Pollinators).
#' betaP,betaA the intraspecific and interspecific competition interaction
#' among species in same groups, Plants group or Animal group.
#' gammaP,gammaA the mutualistic interaction between Plants and Animal.
#' h the Handling time, saturate coefficient.
#' @return the derivation
#' @details .
#' @import deSolve
lv.2 <- function(t, N, parms) {  # Times, State variables, parameters
  alphaP = parms[[1]]
  alphaA = parms[[2]]
  betaP = parms[[3]]
  betaA = parms[[4]]
  gammaP = parms[[5]]
  gammaA = parms[[6]]
  h = parms[[7]]
  
  numP = length(alphaP)
  numA = length(alphaA)
  NP = N[1:numP]
  NAA = N[(numP + 1) : (numP + numA)]
  dNP = alphaP - betaP %*% NP + gammaP %*% NAA / (1 + h * gammaP %*% NAA)
  dNA = alphaA - betaA %*% NAA + gammaA %*% NP / (1 + h * gammaA %*% NP)
  dN = N * c(dNP, dNA)
  list(c(dN))
}

lv.1.check <- function(dataset,  # the network
                       stry = 100,  # number of try for randomizing Structure
                       ptry = 100)  # number of try for Parameterization 
{
  result = data.frame( magnitude = numeric(0), connectance = numeric(0), 
                       randomization = numeric(0), parameterization = numeric(0), nodf = numeric(0),
                       survived = numeric(0), extinct = numeric(0), persistence = numeric(0), lev = numeric(0))
  dataset[dataset > 0] = 1  # get binary network
  numP = dim(dataset)[1]  # number of Plants
  numA = dim(dataset)[2]  # number of Animals
  magnitude = numP * numA  # network magnitude
  connectance = sum(dataset > 0) / magnitude  # network connectance
  
  for (randomization in 1:stry) { 
    dataset.swap = swaplinks(dataset, HowManyToTry = 1000) # randomizing the original network using the swap model
    nodf = nestednodf(dataset.swap)$statistic['NODF']  # nested(dataset.swap, method = 'NODF2')  # NODF of the swapped network
    nest.cmnb = nest.cmnb(dataset.swap)['CMNB']
    for (parameterization in 1:ptry) {
      N0 = runif(numP + numA)  # Initial species abundence, uniformly in (0, 1)
      alphaP = runif(numP) * 0.25 + 0.85  # Intrinsic growth, uniformly in [0.85, 1.1]
      alphaA = runif(numA) * 0.25 + 0.85  # Intrinsic growth, uniformly in [0.85, 1.1]
      
      betaP = matrix(runif(numP * numP), ncol = numP) * 0.02 + 0.22  # Interspecific competition, uniformly in [0.22, 0.24]
      diagP = runif(numP) * 0.02 + 0.99  # Intraspecific competition, uniformly in [0.99. 1.01]
      diag(betaP) = diagP
      betaA = matrix(runif(numA * numA), ncol = numA) * 0.02 + 0.22  # Interspecific competition, uniformly in [0.22, 0.24]
      diagA = runif(numA) * 0.02 + 0.99  # Intraspecific competition, uniformly in [0.99. 1.01]
      diag(betaA) = diagA
      
      gammaP = matrix(runif(numP * numA), ncol = numA) * 0.02 + 0.19  # Mutualistic strength between Plants and Animals, uniformly in [0.19, 0.21]
      gammaP = gammaP * dataset.swap
      gammaA = matrix(runif(numA * numP), ncol = numP) * 0.02 + 0.19  # Mutualistic strength between Plants and Animals, uniformly in [0.19, 0.21]
      gammaA = gammaA * t(dataset.swap)
      
      t = seq(1, 500)  # times
      r = c(alphaP, alphaA)
      M = rbind(cbind(-betaP, gammaP), cbind(gammaA, -betaA))
      parms = list(r, M)  # parameters
      
      lvout <- ode(y = N0, t, func = lv.1, parms = parms)  # , atol = 10^-14, rtol = 10^-12
      lvout.2 <- stode(y = N0, fun = lv.1, parms = parms, pos = TRUE)  # 
      lvout.3 <- runsteady(y = N0, fun = lv.1, parms = parms, times = c(0, 1e5))
      comm = jacobian.full(y = lvout[nrow(lvout), 2:ncol(lvout)], func = lv.1, parms = parms)
      lev = max(Re(eigen(comm)$values))
      
      survived = sum(lvout[nrow(lvout), 2:ncol(lvout)] >= 10^-8)  # Survived species at equillibrim state
      extinct = sum(lvout[nrow(lvout), 2:ncol(lvout)] < 10^-8)  # Extinct species at equillibrim state
      persistence = survived / (numP + numA)  #
      temp = data.frame( magnitude = magnitude, connectance = connectance, 
                         randomization = randomization, parameterization = parameterization, nodf = nodf,
                         survived = survived, extinct = extinct, persistence = persistence, lev = lev)
      result = rbind(result, temp)
    }
  }
  return(result)
}

###############################################################################
### Check the persistence of ODE Equations by Bastolla et al.
lv.2.check <- function(dataset) {
  result = data.frame( magnitude = numeric(0), connectance = numeric(0), 
                       randomization = numeric(0), parameterization = numeric(0), nodf = numeric(0),
                       survived = numeric(0), extinct = numeric(0), persistence = numeric(0), lev = numeric(0))
  dataset[dataset > 0] = 1  # get binary network
  numP = dim(dataset)[1]  # number of Plants
  numA = dim(dataset)[2]  # number of Animals
  magnitude = numP * numA  # network magnitude
  connectance = sum(dataset > 0) / magnitude  # network connectance
  
  for (randomization in 1:250) { 
    dataset.swap = swaplinks(dataset, HowManyToTry = 1000) # randomizing the original network using the swap model
    nodf = nestednodf(dataset.swap)$statistic['NODF']  # nested(dataset.swap, method = 'NODF2')  # NODF of the swapped network
    nest.cmnb = nest.cmnb(dataset.swap)['CMNB']
    for (parameterization in 1:250) {
      N0 = runif(numP + numA)  # Initial species abundence, uniformly in (0, 1)
      alphaP = runif(numP) * 0.25 + 0.85  # Intrinsic growth, uniformly in [0.85, 1.1]
      alphaA = runif(numA) * 0.25 + 0.85  # Intrinsic growth, uniformly in [0.85, 1.1]
      
      betaP = matrix(runif(numP * numP), ncol = numP) * 0.02 + 0.22  # Interspecific competition, uniformly in [0.22, 0.24]
      diagP = runif(numP) * 0.02 + 0.99  # Intraspecific competition, uniformly in [0.99. 1.01]
      diag(betaP) = diagP
      betaA = matrix(runif(numA * numA), ncol = numA) * 0.02 + 0.22  # Interspecific competition, uniformly in [0.22, 0.24]
      diagA = runif(numA) * 0.02 + 0.99  # Intraspecific competition, uniformly in [0.99. 1.01]
      diag(betaA) = diagA
      
      gammaP = matrix(runif(numP * numA), ncol = numA) * 0.02 + 0.19  # Mutualistic strength between Plants and Animals, uniformly in [0.19, 0.21]
      gammaP = gammaP * dataset.swap
      gammaA = matrix(runif(numA * numP), ncol = numP) * 0.02 + 0.19  # Mutualistic strength between Plants and Animals, uniformly in [0.19, 0.21]
      gammaA = gammaA * t(dataset.swap)
      
      h = 0.1  # Saturation parameter
      
      t = seq(1, 500)  # times
      parms = list(alphaP, alphaA, betaP, betaA, gammaP, gammaA, h)  # parameters
      
      lvout <- ode(y = N0, t, func = lv.2, parms = parms)  # , atol = 10^-14, rtol = 10^-12
      lvout.2 <- stode(y = N0, fun = lv.2, parms = parms)  # , pos = TRUE
      lvout.3 <- runsteady(y = N0, fun = lv.2, parms = parms, times = c(0, 1e5))
      comm = jacobian.full(y = lvout[nrow(lvout), 2:ncol(lvout)], func = lv.1, parms = parms)
      lev = max(Re(eigen(comm)$values))
      
      survived = sum(lvout[nrow(lvout), 2:ncol(lvout)] >= 10^-8)  # Survived species at equillibrim state
      extinct = sum(lvout[nrow(lvout), 2:ncol(lvout)] < 10^-8)  # Extinct species at equillibrim state
      persistence = survived / (numP + numA)  #
      temp = data.frame( magnitude = magnitude, connectance = connectance, 
                        randomization = randomization, parameterization = parameterization, nodf = nodf,
                        survived = survived, extinct = extinct, persistence = persistence, lev = lev)
      result = rbind(result, temp)
    }
  }
  return(result)
}

