###############################################################################
# functions accompany with paper:
# Heterogeneity in ecoogical mutualistic networks dominantly determines community stability
# Some codes are extracted and modified from Allesina et al 's codes on:
# https://bitbucket.org/AllesinaLab/superellipses

require('MASS')
library(bipartite)


###############################################################################
#' @title Estimate the maximum eigenvalue of bipartite networks
#' @title using the semicircle (Wigner's law) plus twins method
#' @param A, the incidence matrix of bipartite network
#' @param is.binary, default value is FALSE, that will keep the orginal version of bipartite network,
#'        if is TRUE, will estimate the binary version.
#' @return a list with:
#' @return lev.true (Largest EigenValue), the true largest eigenvalue
#' @return lev.est, the estimated largest eigenvalue
#' @return r, the estimated semicircle radius
#' @return chung, the estimated largest eigenvalue using Chung's method <k^2>/<k>
est.lev.wigner <- function(A, is.binary = FALSE) {
  M = get.adjacency(A)  # get the adjacency matrix of incidence matrix
  if (is.binary)
    M[M != 0] = 1  # convert to binary version
  eigenvalues = eigen(M)$values
  lev.true = max(eigenvalues)  # the true largest eigenvalue
  s = dim(M)[1]  # the number of species (nodes)
  M.tr2 =  sum(diag(M %*% M))  # the trace of second power of adjacency matrix, tr(M^2)
  M.tr4 =  sum(diag(M %*% M %*% M %*% M))  # the trace of fourth power of adjacency matrix, tr(M^4)
  lambda1 = sqrt( ( (s - 2) * M.tr4 + 8 * M.tr2^2 / (s + 2) - 2 * M.tr2^2 ) / (2 * (s + 2)) ) + 
            2 * M.tr2 / (s + 2)
  lambda1 = sqrt(lambda1)
  r = (M.tr2 - 2 * lambda1^2) / (s - 2)
  chung = sum(rowSums(M)**2) / sum(rowSums(M))
  
  # k = sum(M) / s / 2  # Allesina's estimate method for bipartite networks 
  # allesina = sqrt(sqrt(5 * M.tr4 / (s - 1) + 25 * k^2 / s^4) + 5 * k / s^2)
  c(lev.true = lev.true, lev.est = lambda1, r = r, chung = chung)
}

###############################################################################
#' @title Get the adjacency matrix of bipartite network from incidence matrix
#' @param A, the incidence matrix of bipartite network
#' @return Adj, the adjacency matrix of bipartite network
get.adjacency <- function(A){
  NumP <- dim(A)[1]  # number of plants
  NumA <- dim(A)[2]  # number of animals
  S <- NumP + NumA  # number of all species
  Adj <- matrix(0, S, S)  # initialize the adjacency matrix as zero-matrix
  Adj[1:NumP, (NumP + 1):S] <- A  # the upper right sub-matrix is incidence matrix
  Adj <- Adj + t(Adj)  # the lower left sub-matrix is transpose of incidence matrix
  return(Adj)
}


###############################################################################
# Functions of estimating the dominant eigenvalue of bipartite networks 
# using semisuperellipse plus twins method
###############################################################################

###############################################################################
#' @title Estimate the maximum eigenvalue of bipartite networks
#' @title using the semisuperellipse plus twins method
#' @param A, the incidence matrix of bipartite network
#' @param is.binary, default value is FALSE, that will keep the orginal version of bipartite network,
#'        if is TRUE, will estimate the binary version.
#' @return a list with:
#' @return lev.true (Largest EigenValue), the true largest eigenvalue
#' @return lev.est, the estimated largest eigenvalue
#' @return r, the estimated semicircle radius
#' @return chung, the estimated largest eigenvalue using Chung's method <k^2>/<k>
est.lev.semisuperellipse <- function(A, is.binary = FALSE) {
  M = get.adjacency(A)  # get the adjacency matrix of incidence matrix
  if (is.binary)
    M[M != 0] = 1  # convert to binary version
  eigenvalues = eigen(M)$values
  lev.true = max(eigenvalues)  # the true largest eigenvalue
  s = dim(M)[1]  # the number of species (nodes)
  r = eigenvalues[2] / 2
  M.tr2 =  sum(diag(M %*% M))
  M.tr4 =  sum(diag(M %*% M %*% M %*% M))
  M.tr6 =  sum(diag(M %*% M %*% M %*% M %*% M %*% M))
  
  par <- Find.lambda1.and.r.and.n(s, M.tr2, M.tr4, M.tr6, lev.true, r)
  
  c(lev.true = lev.true, lev.est = par[1], r = par[2], n = par[3])
}

#' @title the second raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu2, the second raw moment \mu^2
Mu2 <- function(r,n) {
  mu2 =  r^2 * factorial(2/n) * gamma(3/n) /
    (factorial(1/n) * gamma(4/n))
  return(mu2)
}
Mu2.old <- function(r, n) {  # a more complex version of \mu^2
  mu2 = 4^(1/n) * r^2 * gamma(1/2+1/n) * gamma(3/n) /
    (sqrt(pi) * gamma(4/n))
  return(mu2)
}

#' @title the fourth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu4, the fourth raw moment \mu^4
Mu4 <- function(r,n) {
  mu4 = 8 * r^4 * 2 * gamma(2/n) * gamma(5/n) /
    (3 * gamma(1/n) * gamma(6/n))
  return(mu4)
}

#' @title the sixth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu6, the sixth raw moment \mu^6
Mu6 <- function(r,n) {
  mu4 = 8 * r^6 * 2 * gamma(2/n) * gamma(7/n) /
    ( gamma(1/n) * gamma(8/n))
  return(mu4)
}

#' @title the minimize function for optimization, extracted and modified from Allesina et al 's codes
#' @param pars, the variables be optimized, include \lambda_1, r, n
#' @param s, the number of species
#' @param Tr2, Tr4, Tr6 the trace of second, fourth, sixth power of adjacency matrix respectively
Minimize.for.SuperEllipse <- function(pars, s, Tr2, Tr4, Tr6){
  lambda1 <- pars[1]
  r <- pars[2]
  n <- pars[3]
  mu2 <- Mu2(r,n)
  mu4 <- Mu4(r,n)
  mu6 <- Mu6(r,n)
  ExpectedTr2 <- abs(mu2 * (s-2) + 2*(lambda1)^2)
  ExpectedTr4 <- abs(mu4 * (s-2) + 2*(lambda1)^4)
  ExpectedTr6 <- abs(mu6 * (s-2) + 2*(lambda1)^6)
  term1 <- max(Tr2 / ExpectedTr2, ExpectedTr2 / Tr2)
  term2 <- max(Tr4 / ExpectedTr4, ExpectedTr4 / Tr4)
  term3 <- max(Tr6 / ExpectedTr6, ExpectedTr6 / Tr6)
  return(abs(term1) * abs(term2) * abs(term3))
}

#' @title estimate by optimization, extracted and modified from Allesina et al 's codes
#' @param s, the number of species
#' @param Tr2, Tr4, Tr6 the trace of second, fourth, sixth power of adjacency matrix respectively
#' @param initial guess for \lambda_1, and r
Find.lambda1.and.r.and.n <- function(s, Tr2, Tr4, Tr6, guess.lambda1, guess.r){
  goal <- 100
  guess.n <- seq(0.1, 50, length = 100)
  initial.n <- 1
  best.pars <- NULL
  best.goal <- goal
  ## Use several initial guesses for n, and take the best result
  while (goal > 1.001){
    # print(c(guess.lambda1, guess.r, guess.n[initial.n]))
    z <- optim(c(guess.lambda1, guess.r, guess.n[initial.n]),
               Minimize.for.SuperEllipse,
               s = s,
               Tr2 = Tr2,
               Tr4 = Tr4,
               Tr6 = Tr6
    )
    initial.n <- initial.n + 1
    goal <- z$value
    if (z$par[3] < 0) goal <- 10000
    if (goal < best.goal){
      best.goal <- goal
      best.pars <- z$par
    }
    # print(goal)
    if (initial.n > 100) break
  }
  return(best.pars)
}

# integrate function of superellipse
integrate.func.superellipse <- function(x,r, n, m) {
  x^m * ((2*r)^n - x^n)^(1/n) / 
    (4* (2*r)^2 * gamma(1+1/n)^2 * gamma(1+2/n)^(-1))
}
#integrate function of Wigner's semicircle
integrate.func.wigner <- function(x,r, m) {
  x^m * 2 * ((2*r)^2 - x^2)^(1/2) / 
    ((2*r)^2 * pi)
}

model <- function(params, s, tr2, tr4, tr6) {
  lambda1 = params[1]
  r = params[2]
  n = params[3]
  F1 <- Mu2(r, n) * (s - 2) + 2 * (lambda1)^2 - tr2
  F2 <- Mu4(r, n) * (s - 2) + 2 * (lambda1)^4 - tr4
  F3 <- Mu6(r, n) * (s - 2) + 2 * (lambda1)^6 - tr6
  c(F1, F2, F3)
}

###############################################################################
## Functions for Null Models
## NM1 : Maintain topological structure, Shuffle weights
## NM2 : Maintain topological structure, Reassign weights uniformly
## NM3 : Maintain degree distribution, Shuffle weights
## NM5 : Rewire links uniformly, Shuffle weights
###############################################################################


#' @title generate binormial distribution random variables with fixed sum
#' @param n, number of random variables
#' @param total, the fixed sum
#' @param p, the parameter of binormial distribution
#' @references http://stats.stackexchange.com/questions/14059/generate-uniformly-distributed-weights-that-sum-to-unity
#' http://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
rbinom.fixsum <- function(n, total, p) {
  r = rbinom(n,n,p) # rnorm(n)  # rgamma(n,1,1)
  r = r / sum(r) * total
  return(r)
}

#' @title generate log-normal distribution random variables with fixed sum
#' @param n, number of random variables
#' @param total, the fixed sum
#' @param sigmalog, the standard deviation of log-normal distribution
#' @references #http://www.mathworks.com/matlabcentral/newsreader/view_thread/304141
#' @import MASS
rlnorm.fixsum <- function(n, total, sigmalog) {
  require('MASS')
  mu = total / n
  mulog = log(mu) - sigmalog^2 / 2
  Q = Null(rep(1, n))
  X = (sigmalog * sqrt(n / (n - 1))) * Q %*% rnorm(n - 1);
  X = X + mulog
  return(exp(X))
}

#' @title Get different measures of bipartite networks
#' @param A, the incidence matrix of bipartite network
#' @return a list of:
#' @return 
#' @import bipartite
get.measures <- function(A) {
  # the nestedness measures 
  nestedness = nested(A, method = c('NODF', 'NODF2', 'weighted NODF', 'wine'))
  
  Adj = get.adjacency(A)  # get the adjacency matrix of incidence matrix
  s = dim(Adj)[1]  # the number of species (nodes)
  strengths2 = rowSums(Adj^2)  # square of strength of species
  Hq = sum(strengths2^2)  # the heterogeneity of strength square of species
  Hw = sum(Adj^4)  # the heterogeneity of link weights
  wC4 = sum((Adj %*% Adj)^2) - 2*Hq + Hw  # the number of weighted four-link cycles
  w2 = sum(Adj^2)  #  the heterogeneity of link weights
  w = sum(Adj)  # sum of link weights
  lev.weight = Re(eigen(Adj)$value[1])  # the largest eigenvalue of quantitative bipartite network
  
  Adj[Adj>0] = 1  # get the binary binary
  degrees = rowSums(Adj)  # degree of species
  sk = sum(degrees)  #  edges of network
  Hk = sum(degrees^2)  # the heterogeneity of species degree
  C4 = sum((Adj %*% Adj)^2) - 2*Hk + sk  # the four-link cycles
  lev.bin = Re(eigen(Adj)$value[1])  # the largest eigenvalue of binary bipartite network
  c(s = s, Hq = Hq, Hw = Hw, wC4 = wC4, w2 = w2, w = w, Hk = Hk, C4 = C4, sk = sk, 
    lev.weight = lev.weight, lev.bin = lev.bin, nestedness)
}


###############################################################################
#' @title Swap links of bipartite network, that will keep the node degree distribution.
#'
#' @param B incidence matrix of bipartite network, rows and cols represent two groups of nodes/species
#' @param HowManyToTry the times to try for swapping links
#' @return the incidence matrix whose links being randomly swapped.
#' @details .  
swaplinks <- function(B, HowManyToTry = 5000) {
  count1 <- 0
  NumP <- dim(B)[1]
  NumA <- dim(B)[2]
  while (count1 < HowManyToTry){
    count1 <- count1 + 1
    ## pick two rows and two columns
    row1 <- sample(1:NumP, 1)
    row2 <- sample(1:NumP, 1)
    col1 <- sample(1:NumA, 1)
    col2 <- sample(1:NumA, 1)
    ## check swappable
    if (B[row1, col1] == 0.0 && B[row1, col2] > 0.0 && B[row2, col1] > 0.0 && B[row2, col2] == 0.0){
      ## swap
      B[row1, col1] <- B[row1, col2]
      B[row1, col2] <- 0.0
      B[row2, col2] <- B[row2, col1]
      B[row2, col1] <- 0.0
    }
    else{
      if (B[row1, col1] > 0.0 && B[row1, col2] == 0.0 && B[row2, col1] == 0.0 && B[row2, col2] > 0.0){
        ## swap
        B[row1, col2] <- B[row1, col1]
        B[row1, col1] <- 0.0
        B[row2, col1] <- B[row2, col2]
        B[row2, col2] <- 0.0
      }
    }
  }
  return(B)
}


