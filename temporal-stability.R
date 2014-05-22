library(simecol)

#' @title the parms of mutualistic lv1 model in soft mean field case
#' @param A, the incident matrix of mutualistic networks which are bipartite
#' @param beta0, the mean value of intraspecies competition
#' @param gamma0, the mean value of interspecies cooperation
#' @param alpha0, the intrinsic growh rate
parms.lv1.softmean <- function(A, gamma0 = 1, alpha0 = 1) {
  numP = dim(A)[1]; numA = dim(A)[2]
  beta0 = - max(numP, numA)  # ensure the negative diagonal dominant, and thus feasible fixed point
  require(bipartite)
  M = as.one.mode(A)  # transform to adjacency matrix
  M[M > 0] = gamma0
  diag(M) = beta0
  r = rep(alpha0, numP + numA)
  parmsV = list(r, M)
  initV = - solve(M) %*% r
  list(parmsV, initV)
}

#' @title Lotka-Volterra (LV) Equations of Holling type I
#' @param t time
#' @param N the variables of the LV system, a vector
#' @param parms parameters passed to LV model
#' r the intrinsic growth rate of species, a vector
#' M the interaction matrix between species
#' @return the derivation
#' @details .
#' @import deSolve
model.lv1 <- function(time, init, parms, ...) {
  S = approxTime1(inputs, time, rule=2)["s.in"]
  r = parms[[1]]  # intrinsic growth rate
  M = parms[[2]]  # interaction matrix
  N = init
  dN <- N * (r + M %*% N)
  list(c(dN))
}

parms.and.init = parms.lv1.softmean(Safariland.C.Pre)
LV1 <- odeModel(
  main = model.lv1, 
  parms = parms.and.init[[1]] ,
  times = c(from = 1, to = 50, by = 0.5),
  init = as.numeric(parms.and.init[[2]]),
  solver = 'lsoda')

initfunc(LV1) <- function(obj) {
  tt <- fromtoby(times(obj))
  inputs(obj) <- as.matrix(data.frame(
    time = tt,
    s.in = pmax(rnorm(tt, mean=1, sd=0.5), 0)
  ))
  obj
}



LV1 <- initialize(LV1)
LV1 <- sim(LV1)
#plot(LV1)

