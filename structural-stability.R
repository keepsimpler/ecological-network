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


# Generate the competition matrix M
# First, generate the structure
s = 20
k = 2
G = graph.connected(c(s/2, s/2), k = 2, gtype = 'bipartite')
A = igraph::get.adjacency(G)
A = as.matrix(A)
# Second, appoint values to the links of the structure
gamma0 = -0.2
beta0 = 2
M = A
M[M > 0] = gamma0
diag(M) = beta0
# Third, 
# the eigenvalues and eigenvectors of competition matrix
Lambda = eigen(M)$values
U = eigen(M)$vectors

