source('nestedness.R')
source('DE.R')

n1 = 200  # number of plants
n2 = 200  # number of animal pollinators
k = 2  # average degree of species
G = graph.connected(s = c(n1, n2), k = k, gtype = 'bipartite')  # generate a random connected bipartite graph
A = get.incidence(G)  # get the incidence matrix of bipartite network [G]
numP = dim(A)[1]; numA = dim(A)[2]

#### Check the feasible equilibrium of LV1 model
repeat {  # run ODE untill finding a feasible solution, i.e. all species survived in steady state
  res = lv1.check(A)
  if (res$extinct == 0) break 
}
M = res$parms[[2]]  # interaction matrix of feasible solution
r = res$parms[[1]]  # intrinsic matrix of feasible solution
J = res$comm  # Jacobian matrix at equilibrium of feasible solution
Nstar = res$Nstar  # species abundance at equilibrium
sum(M %*% Nstar + r > 1e-8)  # the relation between species abundance and (M, r) 
sum(Nstar * M - J > 1e-7)  # relation between Jacobian matrix and (M, Nstar)
sum(as.numeric(-solve(M) %*% r) - Nstar > 1e-10)


###############################################################################
# the soft mean field approximation of LV1 model from Bastolla et al. Nature 458, 2009
# check the relation between range of parameters and the stability :
# parameters :
# 1) nestedness : nodf cmnb tr4 
# 2) rtry : number of tried for randomly structure with definite nestedness ( how to keep the nestedness constant?)
# 3) alpha0 : [0.01]
stepwise = 0.05
result = list()
for (i in 1:10) {
  for (alpha0 in seq(stepwise, 1, by = stepwise)) {
    for (beta0 in seq(stepwise, 1, by = stepwise)) {
      for (gamma0 in seq(stepwise, 1, by = stepwise)) {
        alpha0 = 0.3
        beta0 = 0.25
        gamma0 = 0.45
        # N0 = runif(numP + numA) + alpha0  # Initial species abundence, uniformly in [alpha0, alpha0 + 1]
        N0 = c(rep(5, numP + numA - 1), 1)
        res = lv1.check.softmean(dataset = A, alpha0, beta0, gamma0, N0, extinct.threshold = extinct.threshold.default)
        res = c(res, alpha0 = alpha0, beta0 = beta0, gamma0 = gamma0)
        result[[length(result)+1]] = res
        print(paste(alpha0, beta0, gamma0, res$extinct))
      }
    }
  }
}
lvout = res$lvout
t = 1:500
matplot(t, lvout[, 2:51], type = "l", lwd = 1.5)

### Check the hysteresis, two alternative stable states, multiple basins of attractors

df.hysteresis = data.frame(nest.nodf = numeric(0), Nstars = numeric(0), levs = numeric(0))

for (nest in seq(0.1, 0.7, by = 0.05)) {
  ## Generate bipartite network with definate nestedness
  nest = 0.15
  repeat {
    A = rewirelinks.richer(A, 100)
    if (nest.nodf2(A)$NODF > nest) break
  }
  print(nest.nodf2(A)$NODF)
  ## Find a feasible steady state, get its parms
  repeat {
    res = lv2.check(A)
    if (res$extinct == 0) break
  }
  
  parms = res$parms
  Nstar = res$Nstar
  lev = res$lev
  
  ## decrease the intrinsic growth rate of pollinators
  steps = 400
  Nstars = sapply(1:(steps + 1), function(i) rep(0, n1 + n2))
  Nstars[,1] = Nstar
  levs = rep(0, steps + 1)
  levs[1] = lev
  stepwise = 0.01
  alphaA = parms[[2]]
  t = seq(1, 500)  # times
  for (i in 1:steps) {
    parms[[2]] = alphaA - i * stepwise
    lvout <- ode(y = Nstar, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
    Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
    # Nstar[Nstar < 10^-8] = 0  # species with biomass less than the threshold is considered to be extinct
    
    comm = jacobian.full(y = Nstar, func = lv2, parms = parms)  # Jacobian matrix in equilibrium
    lev = max(Re(eigen(comm)$values))
    
    Nstars[, i+1] = Nstar
    levs[i+1] = lev
    print(lev)
  }
  tmp = data.frame(nest.nodf = nest.nodf2(A)$NODF, Nstars = t(Nstars), levs =levs)
  df.hysteresis = rbind(df.hysteresis, tmp)
}

tmp = t(Nstars)
t = seq(1, steps + 1)
matplot(t, tmp[, 1:10], type = 'l')  # plants
matplot(t, tmp[, 11:35], type = 'l')  # animals
lines(t, levs - min(levs), type = 'l')



## decrease the intrinsic growth rate of pollinators
steps = 400
Nstars = sapply(1:(steps * 50), function(i) rep(0, n1 + n2))
levs = rep(0, steps)
stepwise = 0.01
alphaA = parms[[2]]
t = seq(1, 50)  # times
for (i in 1:steps) {
  parms[[2]] = alphaA - i * stepwise
  lvout <- ode(y = Nstar, t, func = lv2, parms = parms, atol = 10^-14, rtol = 10^-12)  # 
  Nstar = lvout[nrow(lvout), 2:ncol(lvout)]  # the species abundance at steady state
  Nstar[Nstar < 10^-8] = 0  # species with biomass less than the threshold is considered to be extinct
  
  comm = jacobian.full(y = Nstar, func = lv2, parms = parms)  # Jacobian matrix in equilibrium
  lev = max(Re(eigen(comm)$values))
  
  Nstars[, ((i-1)*50+1):(i*50)] = lvout[, 2:ncol(lvout)]
  levs[i+1] = lev
  print(lev)
}
tmp = t(Nstars)
t = seq(1, steps * 50)
matplot(t, tmp[, 1:25], type = 'l')  # plants
matplot(t, tmp[, 26:50], type = 'l')  # animals

##### evaluate 'Disentangling nestedness from models of ecological complex'

matplot(t, lvout[, 2:10], type = "l", lwd = 1.5, col = 1)

############ Nestedness definition
source('nestedness.R')  # import the functions defining different nestednesss

A = Safariland.C.Pre  # incident matrix
A[A > 0] = 1  # transform to binary matrix
B = as.one.mode(A)  # transform to adjacency matrix
B2 = B %*% B # get the overlap matrix, which indicate the common neighbors between species
sum(B2)  # the sum of all the common neighbors between species
degrees = rowSums(B)  # the degrees (number of neighbors) of species
degrees^2  # the number of neighbors of neighbors of species
print(paste(sum(B2), sum(degrees^2)))  # 


########## Check the assortativity of random bipartite networks
## 1. the dis-assortativity increase as rewiring links that increase the degrees of node with more degrees (richer get richer)
## 2. swapping links that keep the degree distribution of nodes don't change the assortativity
## 3. when [numP] / [numA] <> 1, the dis-assortativity increase.
numP = 100; numA = 100; k = 3
assort = list()
G = graph.connected(s = c(numP, numA), k = k, gtype = 'bipartite')
for (i in 1:100) {
  print(paste(i, is.connected(G), assortativity.degree(G)))
  A = igraph::get.incidence(G)
  repeat {
    B = rewirelinks.richer(A, HowManyToTry = 50)
    if (is.connected(graph.incidence(B))) {
      G = graph.incidence(B)
      break
    }
  }
  A = igraph::get.incidence(G)
  for (j in 1:10) {
    B2 = swaplinks(A, HowManyToTry = 5000)
    G2 = graph.incidence(B2)
    assort[[(i-1)* 10 + j]] = c(i = (i-1)* 10 + j, 
                                assort1 = assortativity.degree(G), assort2 = assortativity.degree(G2),
                                nest1 = nest.nodf(B)$NODF, nest2 = nest.nodf(B2)$NODF)
  }
}
assort = data.frame(t(as.data.frame.list(assort)))



## degree heterogeneity and assortativity influence on the mean and variance of species abundance at equilibrium.
result.graphs = list()
s1 = s2 = 50
k = 4
G = graph.connected(s = s1, k = k, gtype = 'regular')
A = as.matrix(igraph::get.adjacency(G))  # a regular bipartite random graph
#G2 = graph.incidence(A)
B = A
count = 0
repeat {
  count = count + 1
  shouldcontinue = FALSE
  ## rewiring one link to a random node which has more neighbors
  ## if tring enough times, and still fail to rewire, then [shouldcontinue] is false.
  for (i in 1:5) {
    B = rewirelinks.richer.onestep(B, ntry = 5000)
    if (B$flag == TRUE) {
      shouldcontinue = TRUE
      break
    }
    else {
      B = B$B
    }
  }
  if (!shouldcontinue) break
  B = B$B  # the new graph
  result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = 0, B = B, B2 = B)
  
  ## swapping two links to increase assortativity
  B2 = B
  count.assort = 0
  repeat {
    count.assort = count.assort + 1
    shouldcontinue.assort = FALSE
    for (j in 1:5) {
      B2 = swaplinks.assort.onestep(B2, ntry = 5000)
      if (B2$flag == TRUE) {
        shouldcontinue.assort = TRUE
        break
      }
      else {
        B2 = B2$B
      }
    }
    if (!shouldcontinue.assort) break
    print(paste(count, count.assort))
    B2 = B2$B  # the new graph
    result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = count.assort, B = B, B2 = B2)
  }
  
  ## swapping two links to decrease assortativity
  B2 = B
  count.assort = 0
  repeat {
    count.assort = count.assort - 1
    shouldcontinue.assort = FALSE
    for (j in 1:5) {
      B2 = swaplinks.disassort.onestep(B2, ntry = 5000)
      if (B2$flag == TRUE) {
        shouldcontinue.assort = TRUE
        break
      }
      else {
        B2 = B2$B
      }
    }
    if (!shouldcontinue.assort) break
    print(paste(count, count.assort))
    B2 = B2$B  # the new graph
    result.graphs[[length(result.graphs) + 1]] =  list(count = count, count.assort = count.assort, B = B, B2 = B2)
  }
}

Alpha = runif(s1 + s2, min = 2, max = 2)
beta0 = ceiling( sqrt((s1 + s2) * k) )  # squared root of edges number, to ensure the positive definitive of M
D = diag(rep(beta0, s1 + s2))
tmp = ldply(result.graphs, function(i) {
  B2 = i$B2
  B2 = get.adjacency(B2)
  degrees = rowSums(B2)  # the degrees
  heterogeneity = sum(degrees^2)  # the degree heterogeneity
  degrees2 = rowSums(B2 %*% B2)  # the two-hop degrees
  assortativity = sum( degrees2 / degrees )
  M = D - B2
  Nstar = solve(M) %*% Alpha  # the feasible fixed point
  abundance.mean = mean(Nstar)
  abundance.sd = sd(Nstar)
  c(index = i$count, index.assort = i$count.assort, 
    heterogeneity = heterogeneity, assortativity = assortativity,
    abundance.mean = abundance.mean, abundance.sd = abundance.sd)
})

tmp %.% group_by(index) %.% summarise(length(index.assort))

plot(1:2214, tmp$heterogeneity)
plot(1:2214, tmp$assortativity)
plot(1:2214, tmp$abundance.mean)
plot(1:2214, tmp$abundance.sd)

tmp2 = filter(tmp, index == 40)
tmp2 = select(tmp2, -index)
