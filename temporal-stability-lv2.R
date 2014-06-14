#!/usr/bin/Rscript

library(simecol)
require(plyr)
require(bipartite)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores


#' @title Lotka-Volterra (LV) Equations of Holling type II by Bastolla et al.
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms parameters passed to LV model
#'        r, the intrinsic growth rate of species, a vector
#'        C, the competition matrix in plants and animals
#'        M, the cooperation matrix between plants and animals
#' @param inputs, implicit parameter which intend to replace the [r] 
#'        to reflect the effct of stochastic environments.
#' @return the derivation
#' @details .
#' @import deSolve
model.lv2 <- function(time, init, parms, ...) {
  ## Get data from [inputs] by linearly interpolates according to current [time] step
  S = approxTime1(inputs, time, rule = 2)
  ## Extract the stochastic intrinsic growth rate by removing the [time] column
  r = S[2:length(S)] 
  #r = parms[[1]]  # intrinsic growth rate
  C = parms[[2]]  # the competition matrix
  M = parms[[3]]  # the cooperation matrix
  h = parms[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  list(c(dN))
}

#' @title the [parms] and [init] of mutualistic lv2 model in soft mean field case
#' @param A, the incident matrix of mutualistic networks which are bipartite
#' @param beta0, the mean value of intraspecies competition,
#'        which is endued according to the condition of feasible equilibrium
#' @param gamma0, the mean value of interspecies cooperation
#' @param alpha0, the intrinsic growth rate
#' @param h0, the Handling time, saturate coefficient.
#' @return a list of [parms] and [init] values of [simObj] class
parms.lv2.softmean <- function(A, gamma0 = 1, alpha0 = 1, h0 = 0) {
  numP = dim(A)[1]; numA = dim(A)[2]
  r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  edges = sum(A > 0)  # the number of edges
  # [beta0] take value of squared root of edge number, to ensure the positive definitive of M
  # and thus ensure the only feasible equilibrium.
  beta0 = ceiling( sqrt(edges) )
  C = diag( rep( beta0, numP + numA ) )  # endue the mean value of intraspecies competition
  M = as.one.mode(A)  # transform to adjacency matrix (function of package [bipartite])
  M[M > 0] = gamma0  # endue the mean value of interspecies cooperation
  h = rep(h0, numP + numA)  # endue the mean value of handling time
  parmsV = list(r, C, M, h)  # the [parms] Value of [simObj] class in package [simecol]
  initV = solve(C - M) %*% r  # the [init] Value, which is close to the steady state.
  list(parmsV, initV)
}


#' @title LV2 simulation function
#' @param graphs, 
#' @param alpha0,
#' @param gamma0,
#' @param h0,
#' @param isout, if output the time serials of simulation
#' @param steps and stepwise of simulation
temporal.stability.lv2 <- function(graphs, alpha0 = 1, gamma0 = 1, h0 = 0.01, isout = FALSE, steps = 10000, stepwise = 0.01) {
  LV2 <- odeModel(
    main = model.lv2, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  alpha.mean = alpha0; alpha.deviation = alpha0 / 2  # the mean value and deviation of intrinsic growth rate
  ## the [initnfuc] of LV2 object, 
  ## which is executed during the object creation process that is triggered by [new] or [initialize]  
  initfunc(LV2) <- function(obj) {
    tt <- fromtoby(times(obj))  # time steps
    n <- length(init(obj))  # species number
    t <- length(tt)  # steps number
    #require(plyr)
    # generate the stochastic intrinsic growth rates, and endue to [inputs]
    tmp = ldply(rep(t, n), runif, min = alpha.mean - alpha.deviation, max = alpha.mean + alpha.deviation)
    tmp = data.frame(time = tt, t(tmp))
    inputs(obj) <- as.matrix(tmp)
    obj
  }
  
  result.lv2 = llply(graphs, .parallel = TRUE, function(graph) {
    A = graph$B  # get the incidency matrix
    print(graph$count)
    parms.and.init = parms.lv2.softmean(A, gamma0 = gamma0, alpha0 = alpha0, h0 = h0)
    parms(LV2) = parms.and.init[[1]]
    init(LV2) = as.numeric(parms.and.init[[2]])
    
    LV2 <- sim(LV2)
    
    LV2.out = out(LV2)
    out = LV2.out[2:length(LV2.out)]  
    Sigma = var(out)
    Mean = apply(out, 2, mean)
    if (isout)
      res = list(graph = A, Sigma = Sigma, Mean = Mean, out = LV2.out)
    else
      res = list(graph = A, Sigma = Sigma, Mean = Mean)
    res
  })
  result.lv2
}


# save(result.lv2, file = 'result.lv2.RData')

# result.lv2.stats = ldply(1:100, function(i) {
#   c(i = i, imean.mean = mean(result.lv2[[i]]$Mean), imean.sd = sd(result.lv2[[i]]$Mean), 
#     imean.cv = mean(result.lv2[[i]]$Mean) / sd(result.lv2[[i]]$Mean),
#     icv.mean = mean(result.lv2[[i]]$Mean / sqrt(diag(result.lv2[[i]]$Sigma))))
# })
# 
# plot(result.lv2[[1]]$Mean,  sqrt(diag(result.lv2[[1]]$Sigma)) )
# 
# t = 1:1000
# matplot(t, out[1:1000, 1:100], type = "l", lwd = 0.5)
# 
