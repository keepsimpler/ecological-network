# Bastolla et al. structural stability

# generate a mean-field competition matrix
rho0 = 0.5  
s = 5
M = matrix(rep(rho0, s^2), ncol = s, nrow = s)
diag(M) = 1

# the (rescaled) intrinsic growth rates
Alpha = runif(n = s, min = 0.9, max = 1.1)


# the fixed points
Nstar = solve(M) %*% Alpha
Nstar

# the eigenvalues and eigenvectors of competition matrix
Lambda = eigen(M)$values
U = eigen(M)$vectors

NstarT = t(Nstar) %*% U
AlphaT = t(Alpha) %*% U

AlphaT == (Lambda * NstarT)


# Generate the cooperation matrix M
# First, generate the structure
s1 = 100; s2 = 100
k = 3
G = graph.connected(c(s1, s2), k = k, gtype = 'bipartite')
A = igraph::get.incidence(G)
#A = igraph::get.adjacency(G)

for (i in 1:100) {
  M = as.one.mode(A)
  beta0 = - max(s1, s2) 
  diag(M) = beta0
  U = eigen(M)$vectors
  #Alpha = - U[,1]  # intrinsic growth rate parallel with the principal eigenvector (corresponding to the dominant eigenvalue)
  #Alpha = rep(2, s1 + s2)
  Alpha = runif(s1 + s2, min = 2, max = 2)
  Nstar = - solve(M) %*% Alpha  # the feasible fixed point
  Phi = M * as.vector(Nstar)  # the community matrix
  abund.mean = mean(Nstar)
  abund.sd = sd(Nstar)
  lev = max(eigen(Phi)$values)
  sev = eigen(Phi)$values[s1+s2-1]
  levM = max(eigen(M)$values)
  A = rewirelinks.richer(A, HowManyToTry = 100)
  print(paste(abund.mean, abund.sd, abund.mean / abund.sd, lev, sev, levM))
}


# Second, appoint values to the links of the structure
gamma0 = 1
beta0 = - (s/2 + 0.1)  # ensure the negative diagonal dominant of M, thus the existence of the feasible fixed point
M = A
M[M > 0] = gamma0
diag(M) = beta0
# Third, rewiring links to get more nested structure


Alpha = rep(1, s1 + s2)
Nstar = - solve(M) %*% Alpha
Phi = M * as.vector(Nstar)  # the community matrix
# the eigenvalues and eigenvectors of community matrix
Lambda = eigen(Phi)$values
U = eigen(Phi)$vectors

