library(yuima)
library(plyr)
library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores


grid = setSampling(Terminal = 10, n = 1000)

mod = setModel( drift = "mu - theta * x", diffusion = "sigma", state.var = "x", time.var = "t", solve.var = "x")
X = simulate(mod, xinit = 0, true.parameter = list(theta = 1, sigma = 1), sampling = grid)

simnum = 1000

Xs = ldply(1:simnum, .parallel = TRUE, function(i) {
  X = simulate(mod, xinit = 0, true.parameter = list(theta = 1, sigma = 1, mu = 1), sampling = grid, seed = set.seed(123))
  # c(mean = mean(X@data@original.data), var = var(X@data@original.data))
  t(X@data@original.data)
})

apply(Xs, 2, mean)

## Two dimensional OU process
sol = c('x1', 'x2')
drift = c('-theta11 * x1 - theta12 * x2', '-theta21 * x1 - theta22 * x2')
diffusion = matrix(c('sigma1', '0', '0', 'sigma2'), 2, 2)
mod2 = setModel(drift = drift, diffusion = diffusion, solve.variable = sol)
X2 = simulate(mod2, true.parameter = list(theta11 = 1, theta12 = 0.1, theta21 = 0.1, theta22 = 1, sigma1 = 1, sigma2 = 1),
         sampling = grid)

Xs = ldply(1:simnum, .parallel = TRUE, function(i) {
  #print(i)
  X2 = simulate(mod2, true.parameter = list(theta11 = 1, theta12 = 0.1, theta21 = 0.1, theta22 = 1, sigma1 = 1, sigma2 = 1),
                sampling = grid)
  c(mean = mean(X2@data@original.data), var = var(X2@data@original.data))
})

apply(Xs, 2, mean)
