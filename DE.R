library(deSolve)
library(rootSolve)

extinct.threshold.default = .Machine$double.eps * 100  # threshold of extinction is 100 times of machine precision

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
lv1 <- function(t, N, parms) {
  r = parms[[1]]  # intrinsic growth rate
  M = parms[[2]]  # interaction matrix
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
#' hP,hA the Handling time, saturate coefficient.
#' @return the derivation
#' @details .
#' @import deSolve
lv2 <- function(t, N, parms) {  # Times, State variables, parameters
  alphaP = parms[[1]]
  alphaA = parms[[2]]
  betaP = parms[[3]]
  betaA = parms[[4]]
  gammaP = parms[[5]]
  gammaA = parms[[6]]
  hP = parms[[7]]
  hA = parms[[8]]
  
  numP = length(alphaP)
  numA = length(alphaA)
  NP = N[1:numP]
  NAA = N[(numP + 1) : (numP + numA)]
  dNP = alphaP - betaP %*% NP + gammaP %*% NAA / (1 + hP * gammaP %*% NAA) 
  dNA = alphaA - betaA %*% NAA + gammaA %*% NP / (1 + hA * gammaA %*% NP)
  dN = N * c(dNP, dNA)
  list(c(dN))
}

###############################################################################
#' @title Check the steady state of LV2 model
#' @references 'the sudden collapse of pollinator communities'
#' @param dataset the incidence matrix of bipartite network to be checked
#' @param extinct.threshold used to decide if species has extinct
#' @return a list with:
#' @return survived  the number of survived species when in steady state
#' @return extinct  the number of extinct species when in steady state
#' @return lev  the Largest EiganValue of Jacobian matrix in steady state
#' @details .
lv2.check <- function(dataset, extinct.threshold = 10^-8) {
  dataset[dataset > 0] = 1  # insure input is a binary network
  numP = dim(dataset)[1]  # number of Plants
  numA = dim(dataset)[2]  # number of Animals

  ## Generate Parameters for LV2 model by sampling from uniform distributions
  alphaP = runif(numP) * 0.3 + 0.05  # Intrinsic growth, uniformly in [0.05, 0.35]
  alphaA = runif(numA) * 0.3 + 0.05  # Intrinsic growth, uniformly in [0.05, 0.35]
  
  betaP = matrix(runif(numP * numP), ncol = numP) * 0.04 + 0.01  # Interspecific competition, uniformly in [0.01, 0.05]
  diagP = runif(numP) * 0.3 + 0.8  # Intraspecific competition, uniformly in [0.8, 1.1]
  diag(betaP) = diagP
  betaA = matrix(runif(numA * numA), ncol = numA) * 0.04 + 0.01  # Interspecific competition, uniformly in [0.01, 0.05]
  diagA = runif(numA) * 0.3 + 0.8  # Intraspecific competition, uniformly in [0.8, 1.1]
  diag(betaA) = diagA
  
  gammaP = matrix(runif(numP * numA), ncol = numA) * 0.4 + 0.8  # Mutualistic strength between Plants and Animals, uniformly in [0.8, 1.2]
  gammaP = gammaP * dataset
  gammaA = matrix(runif(numA * numP), ncol = numP) * 0.4 + 0.8  # Mutualistic strength between Plants and Animals, uniformly in [0.8, 1.2]
  gammaA = gammaA * t(dataset)
  
  hP = runif(numP) * 0.15 + 0.15  # Saturation parameter, uniformly in [0.15, 0.3]
  hA = runif(numA) * 0.15 + 0.15  # Saturation parameter, uniformly in [0.15, 0.3]
  
  parms = list(alphaP, alphaA, betaP, betaA, gammaP, gammaA, hP, hA)  # parameters
  t = seq(1, 500)  # times
  N0 = runif(numP + numA) + 0.5  # Initial species abundence, uniformly in [0.5, 1.5]
  
  lvout <- ode(y = N0, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
  Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
  Nstar[Nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
  
  comm = jacobian.full(y = lvout[nrow(lvout), 2:ncol(lvout)], func = lv2, parms = parms)  # Jacobian matrix in equilibrium
  lev = max(Re(eigen(comm)$values))
  
  survived = sum(Nstar > 0)  # Survived species at equillibrim state
  extinct = sum(Nstar == 0)  # Extinct species at equillibrim state
  list(survived = survived, extinct = extinct, lev = lev, parms = parms, Nstar = Nstar)
}

###############################################################################
#' @title Check the steady state of LV1 model
#' @references parameters refered 'the sudden collapse of pollinator communities'
#' @param dataset the incidence matrix of bipartite network to be checked
#' @param extinct.threshold used to decide if species has extinct
#' @return a list with:
#' @return survived  the number of survived species when in steady state
#' @return extinct  the number of extinct species when in steady state
#' @return comm the Jacobian matrix in steady state (two options of Nstar?)
#' @return lev  the Largest EiganValue of Jacobian matrix in steady state
#' @details .
lv1.check <- function(dataset, extinct.threshold = 10^-8) {
  dataset[dataset > 0] = 1  # insure input is a binary network
  numP = dim(dataset)[1]  # number of Plants
  numA = dim(dataset)[2]  # number of Animals

  ## Generate Parameters for LV1 model by sampling from uniform distributions
  alphaP = runif(numP) * 0.3 + 0.05  # Intrinsic growth, uniformly in [0.05, 0.35]
  alphaA = runif(numA) * 0.3 + 0.05  # Intrinsic growth, uniformly in [0.05, 0.35]
  
  betaP = matrix(runif(numP * numP), ncol = numP) * 0.04 + 0.01  # Interspecific competition, uniformly in [0.01, 0.05]
  diagP = runif(numP) * 0.3 + 10.8  # Intraspecific competition, uniformly in [0.8, 1.1]
  diag(betaP) = diagP
  betaA = matrix(runif(numA * numA), ncol = numA) * 0.04 + 0.01  # Interspecific competition, uniformly in [0.01, 0.05]
  diagA = runif(numA) * 0.3 + 10.8  # Intraspecific competition, uniformly in [0.8, 1.1]
  diag(betaA) = diagA
  
  gammaP = matrix(runif(numP * numA), ncol = numA) * 0.4 + 0.8  # Mutualistic strength between Plants and Animals, uniformly in [0.8, 1.2]
  gammaP = gammaP * dataset
  gammaA = matrix(runif(numA * numP), ncol = numP) * 0.4 + 0.8  # Mutualistic strength between Plants and Animals, uniformly in [0.8, 1.2]
  gammaA = gammaA * t(dataset)
  
  r = c(alphaP, alphaA)  # construct the intrinsic growth rate vector
  M = rbind(cbind(-betaP, gammaP), cbind(gammaA, -betaA))  # construct the interaction matrix
  parms = list(r, M)  # parameters
  
  t = seq(1, 500)  # times
  N0 = runif(numP + numA) + 0.5  # Initial species abundence, uniformly in [0.5, 1.5]
  
  lvout <- ode(y = N0, t, func = lv1, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
  Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
  Nstar[Nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
  
  comm = jacobian.full(y = lvout[nrow(lvout), 2:ncol(lvout)], func = lv1, parms = parms)  # Jacobian matrix in equilibrium
  lev = max(Re(eigen(comm)$values))
  
  survived = sum(Nstar > 0)  # Survived species at equillibrim state
  extinct = sum(Nstar == 0)  # Extinct species at equillibrim state
  list(survived = survived, extinct = extinct, lev = lev, parms = parms, Nstar = Nstar, comm = comm)
}

#### the soft mean field approximation of LV1 model from Bastolla et al. Nature 458, 2009
lv1.check.softmean <- function(dataset, alpha0, beta0, gamma0, N0, extinct.threshold = 10^-8) {
  dataset[dataset > 0] = 1  # insure input is a binary network
  numP = dim(dataset)[1]  # number of Plants
  numA = dim(dataset)[2]  # number of Animals
  s = numP + numA  # number of all species
  # the mean value of interspecific interaction strength,
  # the intraspecific interaction strength has been normalized to 1 by \beta_{ij} / \sqrt(\beta_{ii}\beta_{jj})
  #beta0 = 0.2
  betaP = matrix(beta0, ncol = numP, nrow = numP)
  diag(betaP) = 1
  betaA = matrix(beta0, ncol = numA, nrow = numA)
  diag(betaA) = 1
  
  #gamma0 = 0.2  # the mean value of mutualistic interaction strength
  gammaP = gamma0 * dataset
  gammaA = gamma0 * t(dataset)
  
  #alpha0 = 1  # the mean value of intrinsic growth rate
  r = rep(alpha0, numP + numA)  # construct the intrinsic growth rate vector
  M = rbind(cbind(-betaP, gammaP), cbind(gammaA, -betaA))  # construct the interaction matrix
  parms = list(r, M)  # parameters
  
  t = seq(1, 500)  # times
  # N0 = runif(numP + numA) + alpha0  # Initial species abundence, uniformly in [alpha0, alpha0 + 1]
  
  lvout <- ode(y = N0, t, func = lv1, parms = parms, atol = 10^-14, rtol = 10^-12)  # 

  Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
  Nstar[Nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
  
  comm = jacobian.full(y = Nstar, func = lv1, parms = parms)  # Jacobian matrix in equilibrium
  lev = max(Re(eigen(comm)$values))
  
  survived = sum(Nstar > 0)  # Survived species at equillibrim state
  extinct = sum(Nstar == 0)  # Extinct species at equillibrim state
  list(survived = survived, extinct = extinct, lev = lev, parms = parms, Nstar = Nstar, comm = comm, lvout = lvout)
}


###############################################################################
### Check the persistence of ODE Equations by Bastolla et al.
lv2.check.old <- function(dataset,  # the network
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
      N0 = runif(numP + numA) + 0.5  # Initial species abundence, uniformly in [0.5, 1.5]
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

      hP = runif(numP) * 0.1 + 0.05  # Saturation parameter, uniformly in [0.05, 0.15]
      hA = runif(numA) * 0.1 + 0.05  # Saturation parameter, uniformly in [0.05, 0.15]
      
      t = seq(1, 500)  # times
      parms = list(alphaP, alphaA, betaP, betaA, gammaP, gammaA, hP, hA)  # parameters
      
      lvout <- ode(y = N0, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
      extinct.threshold = 10^-8  # 
      Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
      Nstar[Nstar < extinct.threshold] = 0  # species with biomass less than the threshold is considered to be extinct
      
      #lvout.2 <- stode(y = N0, fun = lv2, parms = parms, 
      #                 atol = 10^-14, rtol = 10^-12, pos = TRUE)  # 
      #lvout.3 <- ode(y = lvout.2$y, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
      # lvout.4 <- runsteady(y = N0, fun = lv2, parms = parms, times = c(0, 1e5))
      
      comm = jacobian.full(y = Nstar, func = lv2, parms = parms)
      lev = max(Re(eigen(comm)$values))
      
      survived = sum(Nstar > 0)  # Survived species at equillibrim state
      extinct = sum(Nstar == 0)  # Extinct species at equillibrim state
      persistence = survived / (numP + numA)  #
      temp = data.frame( magnitude = magnitude, connectance = connectance, 
                         randomization = randomization, parameterization = parameterization, nodf = nodf,
                         survived = survived, extinct = extinct, persistence = persistence, lev = lev)
      result = rbind(result, temp)
    }
  }
  return(result)
}

