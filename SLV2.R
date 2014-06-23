#!/usr/bin/Rscript

library(yuima)
library(plyr)
library(bipartite)

library(doMC)  # 
registerDoMC(100)  # register Multi Cores
getDoParWorkers()  # get available Cores

## Stochastic LV2 model

getSolveVariables <- function(m) {
  sol = llply(1:m, function(i) {
    paste('x', i, sep = '')
  })
  as.character(sol)
}

getDrift <- function(m) {
  drift = llply(1:m, function(i) {
    drift1 = paste('x', i, '*(', 'alpha', i, '-', 'theta', i, i, '*x', i, '-', sep = '')
    drift2 = '('
    for (j in 1:m) {
      if (j != i) {
        drift2 = paste(drift2, 'theta', i, j, '*', 'x', j, '+', sep = '')
      }
    }
    substr( drift2, nchar(drift2), nchar(drift2) ) = ')'
    drift2 = paste(drift2, '/(1 - h*', drift2, ')', sep = '')
    drift = paste(drift1, drift2, ')', sep = '')
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

getParameters <- function(m, Alpha, Theta, Sigma, h) {
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
  parameters['h'] = h 
  parameters
}

# return a 3 dimensional array:
# first dimension: simnum, number of simulations
# second dimension: t+1, the simulating steps of one simulation
# third dimension: m, the dimensions of multivariates
sim.slv2 <- function(simnum = 1000, steps = 1000, stepwise = 0.01, m, Alpha, Theta, Sigma, h = 0.1, Xinit) {
  grid = setSampling(Terminal = steps * stepwise, n = steps)
  drift = getDrift(m)
  diffusion = getDiffusion(m)
  solve.variable = getSolveVariables(m)
  mod3 = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, xinit = Xinit)
  parameters = getParameters(m, Alpha, Theta, Sigma, h)
  #X3 = simulate(mod3, true.parameter = parameters, sampling = grid)
  Xs = laply(1:simnum, .parallel = TRUE, function(i) {
    print(i)
    X = simulate(mod3, true.parameter = parameters, sampling = grid)
    matrix(X@data@original.data, ncol = m)
  })
  Xs
}

sim.slv2.graph <- function(graph, alpha0 = 1, beta0 = NULL, gamma0 = 1, sigma0 = 0.05, h = 0.05,
                           simnum = 1000, steps = 1000, stepwise = 0.01) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  if (is.null(beta0)) {
    beta0 = ceiling( sqrt(edges) )
  }
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  A[A > 0] = gamma0
  Theta = D - A  # competition interaction matrix
  Alpha = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  Xinit = solve(Theta) %*% Alpha  # the feasible fixed point
  Sigma = diag(sigma0, numP + numA)
  slv2.out = sim.slv2(simnum = simnum, steps = steps, stepwise = stepwise, m = numP + numA, 
                      Alpha = Alpha, Theta = Theta, Sigma = Sigma, h = h, Xinit = Xinit)  
  slv2.out
}

# m = 3
# Alpha = c(1, 1, 1)
# Theta = matrix(c(1, -0.1, -0.1, -0.1, 1, -0.1, -0.1, -0.1, 1), ncol = m)
# Sigma = diag(0.05, m)
# Xinit = c(1, 1, 1)
# h = 0.01
# slv2.out = sim.slv2(simnum = 1, steps = 1000, m = m, Alpha = Alpha, Theta = Theta, Sigma = Sigma, h = h, Xinit = Xinit)
# matplot(slv2.out, type = 'l', lwd = 0.7)
