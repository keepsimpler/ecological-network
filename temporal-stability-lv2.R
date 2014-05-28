library(simecol)
require(plyr)
require(bipartite)

#' @title Lotka-Volterra (LV) Equations of Holling type II by Bastolla et al.
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms parameters passed to LV model
#' r, the intrinsic growth rate of species, a vector
#' C, the competition matrix in plants and animals
#' M, the cooperation matrix between plants and animals
#' @return the derivation
#' @details .
#' @import deSolve
model.lv2 <- function(time, init, parms, ...) {
  #First, get the vector which linearly interpolates data from [inputs] according to current [time] step
  S = approxTime1(inputs, time, rule = 2)
  r = S[2:length(S)]  # extract the stochastic intrinsic growth rate by removing the [time] column in the vector
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
#' @param beta0, the mean value of intraspecies competition
#' @param gamma0, the mean value of interspecies cooperation
#' @param alpha0, the intrinsic growth rate
#' @param h0, the Handling time, saturate coefficient.
#' @return a list of [parms] and [init] values of [simObj] class
parms.lv2.softmean <- function(A, gamma0 = 1, alpha0 = 1, h0 = 0) {
  numP = dim(A)[1]; numA = dim(A)[2]
  r = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  # ensure negative diagonal dominant matrix, and thus ensure the only and feasible fixed point
  beta0 = max(numP, numA) / 4  
  C = diag(rep(beta0, numP + numA))  # endue the mean value of intraspecies competition
  M = as.one.mode(A)  # transform to adjacency matrix (function of package [bipartite])
  M[M > 0] = gamma0  # endue the mean value of interspecies cooperation
  h = rep(h0, numP + numA)  # endue the mean value of handling time
  parmsV = list(r, C, M, h)  # the [parms] Value of [simObj] class in package [simecol]
  initV = - solve(M - C) %*% r  # the [init] Value, which is close to the steady state.
  list(parmsV, initV)
}

LV2 <- odeModel(
  main = model.lv2, 
  times = c(from = 1, to = 100, by = 0.01),
  solver = 'lsoda')

# the [initnfuc] of LV2 object, which is executed during the object creation process that is triggered by [new] or [initialize]  
alpha.mean = 1; alpha.deviation = 0.2  # the mean value and deviation of intrinsic growth rate
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
result.lv2 = list()
numP = 50; numA = 50; k = 2
G = graph.connected(c(numP, numA), k = k, gtype = "bipartite")
A = igraph::get.incidence(G)

for (i in 1:100) {
  print(i)
  parms.and.init = parms.lv2.softmean(A, h0 = 0.001)
  parms(LV2) = parms.and.init[[1]]
  init(LV2) = as.numeric(parms.and.init[[2]])
  
  LV2 <- sim(LV2)
  
  LV2.out = out(LV2)
  out = LV2.out[2:length(LV2.out)]  
  Sigma = var(out)
  Mean = apply(out, 2, mean)
  result.lv2[[i]] = list(i = i, Sigma = Sigma, Mean = Mean, out = out)
  A = rewirelinks.richer(A, HowManyToTry = 20) 
}

result.lv2.stats = ldply(1:100, function(i) {
  c(i = i, imean.mean = mean(result.lv2[[i]]$Mean), imean.sd = sd(result.lv2[[i]]$Mean), 
    imean.cv = mean(result.lv2[[i]]$Mean) / sd(result.lv2[[i]]$Mean),
    icv.mean = mean(result.lv2[[i]]$Mean / sqrt(diag(result.lv2[[i]]$Sigma))))
})

plot(result.lv2[[1]]$Mean,  sqrt(diag(result.lv2[[1]]$Sigma)) )

t = 1:1000
matplot(t, out[1:1000, 1:100], type = "l", lwd = 0.5)

