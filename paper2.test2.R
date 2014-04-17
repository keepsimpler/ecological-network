source('nestedness.R')
source('DE.R')

n1 = 25  # number of plants
n2 = 25  # number of animal pollinators
k = 4  # average degree of species
G = graph.connected(s = c(n1, n2), k = k, gtype = 'bipartite')  # generate a random connected bipartite graph
A = get.incidence(G)  # get the incidence matrix of bipartite network [G]


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
        res = lv1.check.softmean(dataset = A, alpha0, beta0, gamma0, extinct.threshold = extinct.threshold.default)
        res = c(res, alpha0 = alpha0, beta0 = beta0, gamma0 = gamma0)
        result[[length(result)+1]] = res
        print(paste(alpha0, beta0, gamma0, res$extinct))
      }
    }
  }
}
lvout = res$lvout
matplot(t, lvout[, 2:51], type = "l", lwd = 1.5)

### Check the hysteresis, two alternative stable states, multiple basins of attractors

df.hysteresis = data.frame(nest.nodf = numeric(0), Nstars = numeric(0), levs = numeric(0))

for (nest in seq(0.1, 0.7, by = 0.05)) {
  ## Generate bipartite network with definate nestedness
  #nest = 0.15
  repeat {
    A = rewirelinks.richer(A, 1)
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


