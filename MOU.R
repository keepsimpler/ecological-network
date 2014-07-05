#!/usr/bin/Rscript

library(yuima)
library(plyr)

# library(doMC)  # 
# registerDoMC(12)  # register Multi Cores
# getDoParWorkers()  # get available Cores

## Multivariate OU (Ornstein - Uhlenbeck) process, parameters as matrix

getSolveVariables <- function(m) {
  sol = llply(1:m, function(i) {
    paste('x', i, sep = '')
  })
  as.character(sol)
}

getDrift <- function(m) {
  if (m == 1) return('theta11*x1')
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
sim.mou <- function(simnum = 1000, steps = 1000, stepwise = 0.01, m, Theta, Sigma, Xinit) {
  grid = setSampling(Terminal = steps * stepwise, n = steps)
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

m = 1
Theta = matrix(-1)
Sigma = matrix(0.8)
Xinit = c(0)
mou.out = sim.mou(simnum = 100, steps = 100000, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)

m = 2
Theta = matrix(c(-1, 0.1, 0.8, -1), ncol = m)
Sigma = diag(0.8, m)
Xinit = c(0, 0)
mou.out = sim.mou(simnum = 100, steps = 100000, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
var = diag(Sigma)^2 / 2 * solve(Theta)  # expected variance-covariance matrix for the symmetric [Theta]

m = 3
Theta = matrix(c(-1, 0.1, 0.1, 0.8, -1, 0.1, 0.1, 0.1, -1), ncol = m)
Sigma = diag(0.8, m)
Xinit = c(0, 0, 0)
mou.out = sim.mou(simnum = 100, steps = 100000, stepwise = 0.01, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# matplot(mou.out[2,,], type = 'l')
# mou.XMeans = aaply(mou.out, .margins = c(2, 3), mean)
# mou.XVars = aaply(mou.out, .margins = c(2), var)
# mou.XVars.self = aaply(mou.XVars, .margins = c(1), function(XVar) {
#   c(diag(XVar), XVar[lower.tri(XVar)])
# })
# matplot(mou.XVars.self, type = 'l')
 

# grid = setSampling(Terminal = 10, n = 10000)
# simnum = 100
# 
# ## One dimensional OU process
# mod = setModel( drift = "mu - theta * x", diffusion = "sigma", state.var = "x", time.var = "t", solve.var = "x")
# X = simulate(mod, xinit = 0, true.parameter = list(theta = 1.5, sigma = 0.5, mu = 1), sampling = grid)

# ## Two dimensional OU process
# sol = c('x1', 'x2')
# drift = c('-theta11 * x1 - theta12 * x2', '-theta21 * x1 - theta22 * x2')
# diffusion = matrix(c('sigma1', '0', '0', 'sigma2'), 2, 2)
# mod2 = setModel(drift = drift, diffusion = diffusion, solve.variable = sol)
# X2 = simulate(mod2, true.parameter = list(theta11 = 1, theta12 = 0.1, theta21 = 0.1, theta22 = 1, sigma1 = 1, sigma2 = 1),
#               sampling = grid)
