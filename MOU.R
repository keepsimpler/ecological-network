#!/usr/bin/Rscript

library(yuima)
library(plyr)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores

## Multivariate OU process, parameters as matrix

getSolveVariables <- function(m) {
  sol = llply(1:m, function(i) {
    paste('x', i, sep = '')
  })
  as.character(sol)
}

getDrift <- function(m) {
  drift = llply(1:m, function(i) {
    drift = paste('theta',i,'1*x1', sep = '')
    drift2 = llply(2:m, function(j) {
      paste('+theta',i,j,'*', 'x', j, sep = '')
    })
    drift = paste(drift, paste(drift2, sep = '', collapse = ''), sep = '')
  })
  as.character(drift)
}

getDiffusion <- function(m) {
  diffusion = llply(1:m, function(i) {
    paste('sigma', i, sep = '')  # '*', 'x', i,
  })
  diffusion = as.character(diffusion)
  zeromatrix = matrix(0, ncol = m, nrow = m)
  diag(zeromatrix) = diffusion
  zeromatrix
}

getParameters <- function(m, Theta, Sigma) {
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
  parameters
}

# return a 3 dimensional array:
# first dimension: simnum, number of simulations
# second dimension: t+1, the simulating steps of one simulation
# third dimension: m, the dimensions of multivariates
sim.mou <- function(simnum = 1000, steps = 1000, m, Theta, Sigma, Xinit) {
  grid = setSampling(Terminal = 10, n = steps)
  drift = getDrift(m)
  diffusion = getDiffusion(m)
  solve.variable = getSolveVariables(m)
  mod3 = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, xinit = Xinit)
  parameters = getParameters(m, Theta, Sigma)
  #X3 = simulate(mod3, true.parameter = parameters, sampling = grid)
  Xs = laply(1:simnum, .parallel = TRUE, function(i) {
    print(i)
    X = simulate(mod3, true.parameter = parameters, sampling = grid)
    matrix(X@data@original.data, ncol = m)
  })
  Xs
}

# m = 3
# Theta = - matrix(c(1, 0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), ncol = m)
# Sigma = diag(1, m)
# Xinit = c(1, 1, 1)
# mou.out = sim.mou(simnum = 1, steps = 1000, m = m, Theta = Phi, Sigma = Sigma, Xinit = Xinit)
# matplot(mou.out, type = 'l')



# 
# grid = setSampling(Terminal = 10, n = 1000)
# simnum = 100
# 
# ## One dimensional OU process
# mod = setModel( drift = "mu - theta * x", diffusion = "sigma", state.var = "x", time.var = "t", solve.var = "x")
# X = simulate(mod, xinit = 0, true.parameter = list(theta = 1, sigma = 1, mu = 1), sampling = grid)
# 
# ## Two dimensional OU process
# sol = c('x1', 'x2')
# drift = c('-theta11 * x1 - theta12 * x2', '-theta21 * x1 - theta22 * x2')
# diffusion = matrix(c('sigma1', '0', '0', 'sigma2'), 2, 2)
# mod2 = setModel(drift = drift, diffusion = diffusion, solve.variable = sol)
# X2 = simulate(mod2, true.parameter = list(theta11 = 1, theta12 = 0.1, theta21 = 0.1, theta22 = 1, sigma1 = 1, sigma2 = 1),
#               sampling = grid)
