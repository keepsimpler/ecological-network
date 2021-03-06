#!/usr/bin/Rscript

library(simecol)
require(plyr)
#require(bipartite)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores

#' @title Lotka-Volterra (LV) Equations of Holling type I
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms, parameters passed to LV model:
#'        r, the intrinsic growth rate of species, a vector
#'        M, the interaction matrix between species
#' @param inputs, implicit parameter which intend to replace the [r] 
#'        to reflect the effct of stochastic environments.
#' @return the derivation
#' @details .
#' @import deSolve
model.lv1 <- function(time, init, parms, ...) {
  ## Get data from [inputs] by linearly interpolates according to current [time] step
  S = approxTime1(inputs, time, rule = 2)
  ## Extract the stochastic intrinsic growth rate by removing the [time] column
  r = S[2:length(S)]
  #r = parms[[1]]  # intrinsic growth rate
  M = parms[[2]]  # interaction matrix
  N = init
  dN <- N * (r - M %*% N)
  list(c(dN))
}

#' @title Get the [parms] and [init] of mutualistic lv1 model in soft mean field case
#' @param A, the incident matrix of mutualistic networks which are bipartite
#' @param beta0, the mean value of intraspecies competition,
#'        which is endued according to the condition of feasible equilibrium
#' @param gamma0, the mean value of interspecies cooperation
#' @param alpha0, the intrinsic growth rate
#' @return a list of [parms] and [init] values 
parms.lv1.softmean <- function(A, gamma0 = 1, alpha0 = 1) {
  numP = dim(A)[1]; numA = dim(A)[2]
  edges = sum(A > 0)  # the number of edges
  # [beta0] take value of squared root of edge number, to ensure the positive definitive of M
  # and thus ensure the only feasible equilibrium.
  beta0 = ceiling( sqrt(edges) )
  D = diag( rep( beta0, numP + numA ) )  # endue the mean value of intraspecies competition
  A = inc.to.adj(A)  # transform incidency matrix to adjacency matrix
  A[A > 0] = gamma0  # endue the mean value of interspecies cooperation
  ## the interaction matrix [M] equal to 
  ## the intraspecies competition matrix [D] substract the interspecies cooperation matrix [A]
  M = D - A  
  r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  parmsV = list(r, M)  # the [parms] Value of [simObj] class in package [simecol]
  initV = solve(M) %*% r  # the [init] Value, which is endued the steady state value.
  list(parmsV, initV)
}

#' @title LV1 simulation function
#' @param graphs, 
#' @param alpha0,
#' @param gamma0,
#' @param isout, if output the time serials of simulation
#' @param steps and stepwise of simulation
temporal.stability.lv1 <- function(graphs, alpha0 = 1, gamma0 = 1, isout = FALSE, steps = 10000, stepwise = 0.01) {
  
  LV1 <- odeModel(
    main = model.lv1, 
    times = c(from = 0, to = steps * stepwise, by = stepwise),
    solver = 'lsoda')
  
  alpha.mean = alpha0; alpha.deviation = alpha0 / 2  # the mean value and deviation of intrinsic growth rate
  ## the [initnfuc] of [LV1] object, 
  ## which is executed during the object creation process that is triggered by [new] or [initialize]  
  initfunc(LV1) <- function(obj) {
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
  
  
  result.lv1 = llply(graphs, .parallel = TRUE, function(graph) {
    A = graph$B  # get the incidency matrix
    print(graph$count)
    parms.and.init = parms.lv1.softmean(A, gamma0 = gamma0, alpha0 = alpha0)
    parms(LV1) = parms.and.init[[1]]
    init(LV1) = as.numeric(parms.and.init[[2]])
    
    # LV1 <- initialize(LV1)  # initialize LV1 object which will trigger the execution of [initfunc]
    LV1 <- sim(LV1)
    #plot(LV1)
    
    LV1.out = out(LV1)
    out = LV1.out[2:length(LV1.out)]  
    Sigma = var(out)
    Mean = apply(out, 2, mean)
    if (isout)
      res = list(graph = A, Sigma = Sigma, Mean = Mean, out = LV1.out)
    else
      res = list(graph = A, Sigma = Sigma, Mean = Mean)
    res
  })
  result.lv1
}

#save(result.lv1, file = 'result.lv1.RData')

# for(i in 1:100) {
#   print(i)
#   parms.and.init = parms.lv1.softmean(A)
#   parms(LV1) = parms.and.init[[1]]
#   init(LV1) = as.numeric(parms.and.init[[2]])
#   
#   # LV1 <- initialize(LV1)  # initialize LV1 object which will trigger the execution of [initfunc]
#   LV1 <- sim(LV1)
#   #plot(LV1)
#   
#   LV1.out = out(LV1)
#   out = LV1.out[2:length(LV1.out)]  
#   Sigma = var(out)
#   Mean = apply(out, 2, mean)
#   tsi = apply(out, 2, mean) / apply(out, 2, sd)  # Temporal Stability of Individual species
#   tse = mean(apply(out,1,sum)) / sd(apply(out,1,sum))  # Temporal Stability of Entile community
#   Result2[[i]] = list(i = i, Sigma = Sigma, Mean = Mean, tsi = tsi, tse = tse)
#   A = rewirelinks.richer(A, HowManyToTry = 20)
# }
# 
# tse2 = ldply(1:29, function(i) {
#   Result2[[i]]$tse
# })
# 
# tsi2 = ldply(1:29, function(i) {
#   c(i = i, mean = mean(Result2[[i]]$tsi), sd = sd(Result2[[i]]$tsi), min = min(Result2[[i]]$tsi), max = max(Result2[[i]]$tsi))
#   #Result2[[i]]$tsi
# })
# 
# cv = ldply(1:29, function(i) {
#   c(i = i, cv = mean(Result2[[i]]$Mean) / sd(Result2[[i]]$Mean), mean = mean(Result2[[i]]$Mean), sd = sd(Result2[[i]]$Mean), 
#     var = sum(Result2[[i]]$Sigma), syn = sum(diag(Result2[[i]]$Sigma)) )
# })
# plot(1:100, cv$syn)
# pairs(cv)
# 
# plot(Result2[[1]]$Mean,  sqrt(diag(Result2[[1]]$Sigma)) )
