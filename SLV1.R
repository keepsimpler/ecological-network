#!/usr/bin/Rscript

library(yuima)
library(plyr)
library(bipartite)

#library(doMC)  # 
#registerDoMC(12)  # register Multi Cores
#getDoParWorkers()  # get available Cores

## Stochastic LV1 model

getSolveVariables <- function(m) {
  sol = llply(1:m, function(i) {
    paste('x', i, sep = '')
  })
  as.character(sol)
}

getDrift <- function(m) {
  drift = llply(1:m, function(i) {
    drift1 = paste('alpha', i, '*', 'x', i, sep = '')
    drift2 = llply(1:m, function(j) {
      paste('-theta',i,j,'*', 'x', i, '*', 'x', j, sep = '')
    })
    drift = paste(drift1, paste(drift2, sep = '', collapse = ''), sep = '')
  })
  as.character(drift)
}

getDiffusion <- function(m) {
  diffusion = llply(1:m, function(i) {
    paste('sigma', i, '*', 'x', i, sep = '')
  })
  diffusion = as.character(diffusion)
  zeromatrix = matrix('0', ncol = m, nrow = m)
  diag(zeromatrix) = diffusion
  zeromatrix
}

getParameters <- function(m, Alpha, Theta, Sigma) {
  parameters = list()
  for(i in 1:m) {
    for(j in 1:m) {
      paraname = paste('theta', i, j, sep = '')
      parameters[paraname] = Theta[i, j]
    }
  }
  for(i in 1:m) {
    paraname = paste('sigma', i, sep = '')
    parameters[paraname] = Sigma[i, i]
  }
  for(i in 1:m) {
    paraname = paste('alpha', i, sep = '')
    parameters[paraname] = Alpha[i]
  }
  parameters
}

# return a 3 dimensional array:
# first dimension: simnum, number of simulations
# second dimension: t+1, the simulating steps of one simulation
# third dimension: m, the dimensions of multivariates
sim.slv1 <- function(simnum = 1000, steps = 1000, stepwise = 0.01, m, Alpha, Theta, Sigma, Xinit) {
  grid = setSampling(Terminal = steps * stepwise, n = steps)
  drift = getDrift(m)
  diffusion = getDiffusion(m)
  solve.variable = getSolveVariables(m)
  mod3 = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, xinit = Xinit)
  parameters = getParameters(m, Alpha, Theta, Sigma)
  #X3 = simulate(mod3, true.parameter = parameters, sampling = grid)
  Xs = laply(1:simnum, .parallel = TRUE, function(i) {
    print(i)
    X = simulate(mod3, true.parameter = parameters, sampling = grid)
    matrix(X@data@original.data, ncol = m)
  })
  Xs
}

sim.slv1.graph <- function(graph, alpha0 = 1, beta0 = 1, gamma0 = 0.1, sigma0 = 0.05,
                           simnum = 1000, steps = 1000, stepwise = 0.01) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
#  if (is.null(beta0)) {
#    beta0 = ceiling( sqrt(edges) )
#  }
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  A[A > 0] = gamma0
  Theta = D - A  # competition interaction matrix
  Alpha = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  Xinit = solve(Theta) %*% Alpha  # the feasible fixed point
  Sigma = diag(sigma0, numP + numA)
  slv1.out = sim.slv1(simnum = simnum, steps = steps, stepwise = stepwise, m = numP + numA, 
                      Alpha = Alpha, Theta = Theta, Sigma = Sigma, Xinit = Xinit)  
  slv1.out
}



## Analysis results for the Stochastic LV model
analysis.slv1.graph <- function(graph, alpha0 = 1, beta0 = 1, gamma0 = NULL, sigma0 = 0.1) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  if (is.null(gamma0)) {
    gamma0 = 1 / sqrt(edges)
  }
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  A[A > 0] = gamma0
  Theta = D - A  # competition interaction matrix
  Alpha = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  Xinit = solve(Theta) %*% Alpha  # the feasible fixed point
  Sigma = diag(sigma0, numP + numA)
  means = (Alpha - diag(Sigma)^2/2) %*% solve(Theta)  # expected means
  vars = c(mu) * diag(Sigma)^2/2 * solve(Theta)
  list(means = means, vars = vars)
}

## get the max cooperation strength [gamma0], that allow a feasible steady state exists
get.gamma0.max <- function(graph, beta0 = 1) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  gamma0 = 1 / sqrt(edges) 
  repeat {
    gamma0 = gamma0 + 0.0001
    A[A > 0] = gamma0
    Theta = D - A  # competition interaction matrix
    if (min(eigen(Theta)$values) <= 0.01) {
      gamma0 = gamma0 - 0.0001
      break
    }
  }
  gamma0
}




# m = 1
# Alpha = c(1)
# Theta = matrix(1)
# Sigma = matrix(0.8)
# Xinit = c(1)
# slv1.out = sim.slv1(simnum = 100, steps = 100000, m = m, Alpha = Alpha, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# 
# m = 2
# Alpha = c(1, 1)
# Theta = matrix(c(1, -0.1, -0.8, 1), ncol = m)
# Sigma = diag(0.8, m)
# Xinit = c(1, 1)
# slv1.out = sim.slv1(simnum = 100, steps = 100000, m = m, Alpha = Alpha, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# mu = (Alpha - diag(Sigma)^2/2) %*% solve(Theta)  # expected mean
# 
# m = 3
# Alpha = c(1, 1, 1)
# Theta = matrix(c(1, -0.1, -0.1, -0.8, 1, -0.1, -0.1, -0.1, 1), ncol = m)
# Sigma = diag(0.8, m)
# Xinit = c(1, 1, 1)
# slv1.out = sim.slv1(simnum = 100, steps = 100000, m = m, Alpha = Alpha, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# matplot(slv1.out, type = 'l', lwd = 0.7)
