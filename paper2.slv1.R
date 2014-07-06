library(plyr)

slv1 <- function(steps = 1000, stepwise = 0.01, Alpha, Theta, sigma0, Xinit, m) {
  dW = mvrnorm(steps, mu = rep(0, m), Sigma = diag(sqrt(stepwise), m))
  X = array(dim = c(steps + 1, m))
  X[1, ] = Xinit
  for (i in 2:(steps + 1)) {
    X[i, ] = X[i - 1] + X[i - 1] * (Alpha - Theta %*% X[i - 1, ] + sigma0 * dW[i - 1, ])
  }
  X
}

slv1.graph <- function(graph, alpha0 = 1, beta0 = 1, gamma0 = 0.1, sigma0 = 0.05,
                           simnum = 100, steps = 10000, stepwise = 0.01) {
  edges = sum(graph > 0)
  numP = dim(graph)[1]
  numA = dim(graph)[2]
  D = diag(beta0, numP + numA)
  A = as.one.mode(graph)
  A[A > 0] = gamma0
  Theta = D - A  # competition interaction matrix
  Alpha = rep(alpha0, numP + numA)  # endue the mean value of intrinsic growth rate
  Xinit = solve(Theta) %*% Alpha  # the feasible fixed point
  Sigma = diag(sigma0, numP + numA)
  slv1.out = laply(1:simnum, function(i) {
    print(i)
    slv1(steps = steps, stepwise = stepwise, Alpha = Alpha, Theta = Theta, 
         sigma0 = sigma0, Xinit = Xinit, m = numP + numA)
  })
  slv1.out
}

