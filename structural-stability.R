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

library(bipartite)

# Generate the cooperation matrix M
# First, generate the structure
s1 = 100; s2 = 100
k = 3
G = graph.connected(c(s1, s2), k = k, gtype = 'bipartite')
A = igraph::get.incidence(G)
#A = igraph::get.adjacency(G)

#ss = data.frame(i = numeric(0), abund.mean = numeric(0), abund.sd = numeric(0), lev = numeric(0), sev = numeric(0),
#                nlev = numeric(0), nsev = numeric(0), eigs.prod = numeric(0), degree.mean = numeric(0), degree.sd = numeric(0))  # Structural Stability
ss = list()
for (i in 1:150) {
  M = as.one.mode(A)
  M2 = M %*% M
  degrees = rowSums(M)
  degrees2 = rowSums(M2)
  #degree.mean = mean(degrees)
  #degree.sd = sd(degrees)
  beta0 = - max(s1, s2)/4
  #beta0 = - rowSums(M) + 1
  diag(M) = beta0
  U = eigen(M)$vectors
  eigs = eigen(M)$values
  Alpha = runif(s1 + s2, min = 2, max = 2)
  Nstar = - solve(M) %*% Alpha  # the feasible fixed point
  NstarT = t(Nstar) %*% U
  AlphaT = t(Alpha) %*% U
  AlphaT == (eigs * NstarT)
  #Alpha = - U[,1]  # intrinsic growth rate parallel with the principal eigenvector (corresponding to the dominant eigenvalue)
  #Alpha = rep(2, s1 + s2)
  #Phi = M * as.vector(Nstar)  # the community matrix
  #abund.mean = mean(Nstar)
  #abund.sd = sd(Nstar)
  #lev = eigs[1]
  #sev = eigs[2]
  #nlev = eigs[s1 + s2]
  #nsev = eigs[s1 + s2 - 1]
  #eigs.sum = sum(eigs)
  #eigs.prod = prod(eigs)^(1/(s1+s2))
  #levM = max(eigen(M)$values)
  # Third, rewiring links to get more nested structure
  A = rewirelinks.richer(A, HowManyToTry = 50)
  #tmp = data.frame(i = i, abund.mean = abund.mean, abund.sd = abund.sd, lev = lev, sev = sev, nlev = nlev, 
  #                 nsev = nsev, eigs.prod = eigs.prod, degree.mean = degree.mean, degree.sd = degree.sd)
  #ss = rbind(ss, tmp)
  ss[[i]] = list(eigs = eigs, Nstar = Nstar, degrees = degrees, degrees2 = degrees2, 
                 U = U, NstarT = NstarT, AlphaT = AlphaT)
  #print(paste(abund.mean, abund.sd, abund.mean / abund.sd, lev, sev, nlev, nsev, eigs.sum, eigs.prod))
}

tmp = ldply(ss, function(i) {
  c(Nstar.mean = mean(i$Nstar), Nstar.sd = sd(i$Nstar), Nstar.shannon = sum(i$Nstar * log(i$Nstar)), eigs.sd = sd(i$eigs), 
    eiv.mean = mean(i$AlphaT * colSums(i$U) / i$eigs), eiv.sd = sd(i$AlphaT * colSums(i$U) / i$eigs),
    eiv.max = min(colSums(i$U)), eig.max = min(i$eigs), eig.sev =i$eigs[200]/i$eigs[199])
})

ss[[1]]$Nstar[1] == - sum(ss[[1]]$AlphaT * ss[[1]]$U[1,] / ss[[1]]$eigs)
sum(colSums(ss[[1]]$U) * ss[[1]]$AlphaT) / (s1 + s2) == mean(Alpha)
sum(ss[[1]]$AlphaT * colSums(ss[[1]]$U) / ss[[1]]$eigs)

plot(ss[[50]]$AlphaT * colSums(ss[[50]]$U), ss[[50]]$eigs )

colSums(ss[[1]]$U) / ss[[1]]$AlphaT = Alpha


# check the relation between eigvectors and eigenvalues
s = 200; k = 3
G = graph.connected(s, k, gtype = 'er')
A = igraph::get.adjacency(G)
A = as.matrix(A)
#D = diag(rowSums(A) + 0.1) 
D = diag(rep(s / 10, s))
L = D - A
#Lambda = eigen(L)$values
U = eigen(L)$vectors
eigs = eigen(L)$values
Alpha = runif(s, min = 2, max = 2)
Nstar = solve(L) %*% Alpha  # the feasible fixed point
NstarT = t(Nstar) %*% U
AlphaT = t(Alpha) %*% U
AlphaT == (eigs * NstarT)

AlphaT[200] * U[,200] / eigs[200]



# check the relation between distribution of species abundances and the degrees and second-order degrees
s1 = 100; s2 = 100
k = 3
G = graph.connected(c(s1, s2), k = k, gtype = 'bipartite')
A = igraph::get.incidence(G)

V = - as.one.mode(A)
beta0 = ceiling( sqrt((s1 + s2) * k) )  # squared root of edges number, to ensure the positive definitive of M
D = diag(rep(beta0, s1 + s2))
M = D + V
Drev = solve(D)
I = diag(rep(1, s1 + s2))
(I + V %*% Drev) %*% D == M
B = V %*% Drev
Mrev = Drev - Drev^2 %*% V + Drev^3 %*% V %*% V - Drev^4 %*% V %*% V %*% V + Drev^5 %*% V %*% V %*% V %*% V
(solve(M) - Mrev)

