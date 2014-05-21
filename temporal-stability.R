# The temporal stability of Multivariate AR(1) model


#' @title Multivariate AR(1) model simulation
#' @references IVES, ESTIMATING COMMUNITY STABILITY AND ECOLOGICAL INTERACTIONS FROM TIME-SERIES DATA
#' @references X[t] = A + B * X[t - 1] + E[t]
#' @param A, the intrinsic growh rate
#' @param B, the interaction matrix
#' @param SigmaE, the covariance matrix of noise (errors)
#' @param t, times
#' @param X0, the initial state
mar1 <- function(A, B, SigmaE, t, X0 = NULL) {
  if (dim(B)[1] != dim(B)[2])
    stop('Interaction matrix B should be a square matrix.')
  m = dim(B)[1]  # dimension of AR model
  eig = eigen(B, only.values = TRUE)$values
  if ((any(Mod(eig) > 1)))
    warning('Unstable AR model.')
  require(MASS)
  # m*1 vector of process errors that has a multivariate normal distribution
  # with mean vector 0 and covariance matrix SigmaE
  E = mvrnorm(t, rep(0, m), SigmaE)  
  E = rep(1, t) %*% t(A) + E  # A + E
  I = diag(1, m)
  Mu = solve(I - B) %*% A  # mean vector of stationary distribution
  Mu = as.vector(Mu)
  if (is.null(X0))  # if initial state values are not provided, use the mean vector of stationary distribution.
    X0 = Mu
  X = matrix(0, nrow = 1 + t, ncol = m)
  X[1, ] = X0
  for (k in seq(1, t))
    X[k + 1, ] = X[k, ] %*% B + E[k, ]
  X
}

A=c(1,1)
B=rbind(c(0.7, -0.7),c(-0.7,0.7))
m = dim(B)[1]
SigmaE=rbind(c(1,0),c(0,1))
t = 10000
out = mar1(A, B, SigmaE, t, X0 = c(0,0))
apply(out[1000:t,], 2, mean)
var(out[1000:t,])
solve(diag(1,m) - B) %*% A
solve(diag(1,m^2) - kronecker(B, B)) %*% as.vector(SigmaE)
plot((t+1-1000):(t+1), out[(t+1-1000):(t+1),2], type = 'l')

sum(apply(out,2,mean)) == mean(apply(out,1,sum)) 
sum(var(out)) == var(apply(out,1,sum))

MeanX = mean(apply(out,1,sum)) 
tse = mean(apply(out,1,sum)) / sd(apply(out,1,sum))  # Temporal Stability of Entile community
tsi = mean(apply(out, 2, mean) / apply(out, 2, sd))  # Temporal Stability of Individual species
SigmaX = var(out)
VarX = sum(diag(SigmaX))
CovX = sum(SigmaX) - VarX

v = -0.5
m = 5
B = matrix(rep(v, m^2), ncol = m, nrow = m)
diag(B) = 1
I = diag(1, m)
solve(I - B)

library(plyr)
ldply(1:15, function(m) {
  v = -0.5
  B = matrix(rep(v, m^2), ncol = m, nrow = m)
  I = diag(1, m)
  A = (I - B) %*% rep(1, m)
  SigmaE = diag(1, m)
  t = 2000
  out = mar1(A, B, SigmaE, t)
  out  
})

for (m in 1:15) {
  v = -0.5
  B = matrix(rep(v, m^2), ncol = m, nrow = m)
  A = (I - B) %*% rep(1, m)
  SigmaE = diag(1, m)
  t = 2000
  out = mar1(A, B, SigmaE, t)
  
}



