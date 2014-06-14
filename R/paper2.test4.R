graph = result.graphs.0[[126]]$B  # 
edges = sum(graph > 0)
numP = dim(graph)[1]
numA = dim(graph)[2]
beta0 = ceiling( sqrt(edges) )
beta0 = max(numP, numA)
D = diag(beta0, numP + numA)
A = as.one.mode(graph)
M = D - A  # competition interaction matrix
Alpha = runif(numP + numA, min = 1, max = 1)
Nstar = solve(M) %*% Alpha  # the feasible fixed point
Phi = - M * as.vector(Nstar)  # the community matrix

## AR(1) test
A = c(1)
B = matrix(0.8)
m = dim(B)[1]
SigmaE = matrix(1)
t = 10000
out = mar1(A, B, SigmaE, t, X0 = c(1))
