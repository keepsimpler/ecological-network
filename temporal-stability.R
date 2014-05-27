library(simecol)
require(plyr)

library(simecol)

#' @title the [parms] and [init] of mutualistic lv1 model in soft mean field case
#' @param A, the incident matrix of mutualistic networks which are bipartite
#' @param beta0, the mean value of intraspecies competition
#' @param gamma0, the mean value of interspecies cooperation
#' @param alpha0, the intrinsic growth rate
#' @return a list of [parms] and [init] values of [simObj] class
parms.lv1.softmean <- function(A, gamma0 = 1, alpha0 = 1) {
  numP = dim(A)[1]; numA = dim(A)[2]
  beta0 = - max(numP, numA) / 4  # ensure negative diagonal dominant matrix, and thus ensure the only and feasible fixed point
  require(bipartite)
  M = as.one.mode(A)  # transform to adjacency matrix
  M[M > 0] = gamma0  # endue the mean value of interspecies cooperation
  diag(M) = beta0  # endue the mean value of intraspecies competition
  r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  parmsV = list(r, M)  # the [parms] Value of [simObj] class in package [simecol]
  initV = - solve(M) %*% r  # the [init] Value
  list(parmsV, initV)
}

#' @title Lotka-Volterra (LV) Equations of Holling type I
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms parameters passed to LV model
#' r the intrinsic growth rate of species, a vector
#' M the interaction matrix between species
#' @return the derivation
#' @details .
#' @import deSolve
model.lv1 <- function(time, init, parms, ...) {
  #First, get the vector which linearly interpolates data from [inputs] according to current [time] step
  S = approxTime1(inputs, time, rule=2)
  r = S[2:length(S)]  # extract the stochastic intrinsic growth rate by removing the [time] column in the vector
  #r = parms[[1]]  # intrinsic growth rate
  M = parms[[2]]  # interaction matrix
  N = init
  dN <- N * (r + M %*% N)
  list(c(dN))
}


# Result = llply(5:15, function(i) {
#   print(i)
#   numP = i; numA = i; k = 2
#   G = graph.connected(c(numP, numA), k = k, gtype = "bipartite")
#   A = igraph::get.incidence(G)
#   
#   parms.and.init = parms.lv1.softmean(A)
#   LV1 <- odeModel(
#     main = model.lv1, 
#     parms = parms.and.init[[1]] ,
#     times = c(from = 1, to = 50, by = 0.01),
#     init = as.numeric(parms.and.init[[2]]),
#     solver = 'lsoda')
#   
#   # the [initnfuc] of LV1 object, which is executed during the object creation process that is triggered by [new] or [initialize]  
#   alpha.mean = 1; alpha.deviation = 0.5  # the mean value and deviation of intrinsic growth rate
#   initfunc(LV1) <- function(obj) {
#     tt <- fromtoby(times(obj))  # time steps
#     n <- length(init(obj))  # species number
#     t <- length(tt)  # steps number
#     require(plyr)
#     # generate the stochastic intrinsic growth rates, and endue to [inputs]
#     tmp = ldply(rep(n, t), runif, min = alpha.mean - alpha.deviation, max = alpha.mean + alpha.deviation)
#     tmp = data.frame(time = tt, tmp)
#     inputs(obj) <- as.matrix(tmp)
#     obj
#   }
#   
#   
#   LV1 <- initialize(LV1)  # initialize LV1 object which will trigger the execution of [initfunc]
#   LV1 <- sim(LV1)
#   #plot(LV1)
#   
#   LV1.out = out(LV1)
#   out = LV1.out[2:length(LV1.out)]  
#   Sigma = var(out)
#   Mean = apply(out, 2, mean)
#   tsi = apply(out, 2, mean) / apply(out, 2, sd)  # Temporal Stability of Individual species
#   tse = mean(apply(out,1,sum)) / sd(apply(out,1,sum))  # Temporal Stability of Entile community
#   list(i = i, Sigma = Sigma, Mean = Mean, tsi = tsi, tse = tse)
# })
# 
# tse = ldply(1:11, function(i) {
#   Result[[i]]$tse
# })
# 
# tsi = ldply(1:11, function(i) {
#   c(mean = mean(Result[[i]]$tsi), sd = sd(Result[[i]]$tsi))
# })
# 
# plot(seq(1,11), tsi$sd)
# 


LV1 <- odeModel(
  main = model.lv1, 
  times = c(from = 1, to = 100, by = 0.01),
  solver = 'lsoda')

# the [initnfuc] of LV1 object, which is executed during the object creation process that is triggered by [new] or [initialize]  
alpha.mean = 1; alpha.deviation = 0.5  # the mean value and deviation of intrinsic growth rate
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


Result2 = list()
numP = 50; numA = 50; k = 2
G = graph.connected(c(numP, numA), k = k, gtype = "bipartite")
A = igraph::get.incidence(G)
for(i in 1:100) {
  print(i)
  parms.and.init = parms.lv1.softmean(A)
  parms(LV1) = parms.and.init[[1]]
  init(LV1) = as.numeric(parms.and.init[[2]])
  
  # LV1 <- initialize(LV1)  # initialize LV1 object which will trigger the execution of [initfunc]
  LV1 <- sim(LV1)
  #plot(LV1)
  
  LV1.out = out(LV1)
  out = LV1.out[2:length(LV1.out)]  
  Sigma = var(out)
  Mean = apply(out, 2, mean)
  tsi = apply(out, 2, mean) / apply(out, 2, sd)  # Temporal Stability of Individual species
  tse = mean(apply(out,1,sum)) / sd(apply(out,1,sum))  # Temporal Stability of Entile community
  Result2[[i]] = list(i = i, Sigma = Sigma, Mean = Mean, tsi = tsi, tse = tse)
  A = rewirelinks.richer(A, HowManyToTry = 20)
}

tse2 = ldply(1:29, function(i) {
  Result2[[i]]$tse
})

tsi2 = ldply(1:29, function(i) {
  c(i = i, mean = mean(Result2[[i]]$tsi), sd = sd(Result2[[i]]$tsi), min = min(Result2[[i]]$tsi), max = max(Result2[[i]]$tsi))
  #Result2[[i]]$tsi
})

cv = ldply(1:29, function(i) {
  c(i = i, cv = mean(Result2[[i]]$Mean) / sd(Result2[[i]]$Mean), mean = mean(Result2[[i]]$Mean), sd = sd(Result2[[i]]$Mean), 
    var = sum(Result2[[i]]$Sigma), syn = sum(diag(Result2[[i]]$Sigma)) )
})
plot(1:100, cv$syn)
pairs(cv)

plot(Result2[[1]]$Mean,  sqrt(diag(Result2[[1]]$Sigma)) )
