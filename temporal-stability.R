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
  list(r, M)
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
model.lv1 <- function(t, N, parms) {
  r = parms[[1]]  # intrinsic growth rate
  M = parms[[2]]  # interaction matrix
  dN <- N * (r + M %*% N)
  list(c(dN))
}


LV1 <- odeModel(
  main = model.lv1, 
  parms = parms.lv1.softmean(Safariland.C.Pre),
  times = c(from=1, to=50, by = 1),
  init = rep(1, 31),
  solver = 'lsoda')

LV1 <- sim(LV1)
plot(LV1)

