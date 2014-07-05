################################################################################
# Analyze and display the results

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

## Input: A four dimensional array, which is the output of simulation of Stochastic LV2 model
## First dimension: the index of graphs with increasing nestedness
## Second dimension: the number of simulation of Stochastic LV2 model
## Third dimension: the time steps of simulation
## Fourth dimension: the number of species

## Get the mean vector of Stochastic LV2 model from the output
#slv1.XMeans = aaply(slv1.out, .margins = c(1,3,4), mean) # three dimensional matrix: graphs * (t+1) * m
slv1.XMeans = array(dim = c(2, 1001, 50))
slv1.XMeans[1,,] = slv1.1.XMeans
slv1.XMeans[2,,] = slv1.124.XMeans
slv1.XMeans = slv1.XMeans[, 502:1001,] # get the subset of stationary distribution
slv1.XMeans.mean = aaply(slv1.XMeans, .margins = c(1, 3), mean)
slv1.XMeans.mean.mean = aaply(slv1.XMeans.mean, 1, mean)
slv1.XMeans.mean = data.frame(slv1.XMeans.mean)
slv1.XMeans.mean$i = 1:dim(slv1.XMeans)[1]
slv1.XMeans.mean.long = melt(slv1.XMeans.mean, id.vars = 'i')
slv1.XMeans.mean.mean = data.frame(slv1.XMeans.mean.mean)
slv1.XMeans.mean.mean$i = 1:dim(slv1.XMeans)[1]
slv1.XMeans.mean.mean.long = melt(slv1.XMeans.mean.mean, id.vars = 'i')
p <- ggplot(slv1.XMeans.mean.long, aes(x = factor(i), y = value, colour = factor(i))) +
  geom_point(alpha = 1, size = 1) + 
  geom_point(data = slv1.XMeans.mean.mean.long, aes(x = factor(i), y = value), alpha = 1, size = 3) +  #stat = 'summary', fun.y = 'mean',
  scale_y_continuous( breaks = c(seq(1, 3, by = 0.5), 4, 5, 10) ) +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Biomass of species') +
  theme_bw()
p


## Get the variance-covariance matrix \Sigma of stochastic LV2 model from the output
#slv1.XVars = aaply(slv1.out, .margins = c(1,3), var) # four dimensional matrix: graphs * (t+1) * m * m
slv1.XVars = array(dim = c(2, 1001, 50, 50))
slv1.XVars[1,,,] = slv1.1.XVars
slv1.XVars[2,,,] = slv1.124.XVars
slv1.XVars = slv1.XVars[, 502:1001,,] # get the subset of stationary distribution

## the variances of individual species
slv1.XVars.self = aaply(slv1.XVars, .margins = c(1, 2), function(XVar) {
  c(diag(XVar))
})
slv1.XVars.self.mean = aaply(slv1.XVars.self, .margins = c(1, 3), mean)
slv1.XVars.self.mean.mean = aaply(slv1.XVars.self.mean, 1, mean)
slv1.XVars.self.mean = data.frame(slv1.XVars.self.mean)
slv1.XVars.self.mean$i = 1:dim(slv1.XVars)[1]
slv1.XVars.self.mean.long = melt(slv1.XVars.self.mean, id.vars = 'i')
slv1.XVars.self.mean.mean = data.frame(slv1.XVars.self.mean.mean)
slv1.XVars.self.mean.mean$i = 1:dim(slv1.XVars)[1]
slv1.XVars.self.mean.mean.long = melt(slv1.XVars.self.mean.mean, id.vars = 'i')
p <- ggplot(slv1.XVars.self.mean.long, aes(x = factor(i), y = value, colour = factor(i))) +
  geom_point(alpha = 1, size = 1) + 
  geom_point(data = slv1.XVars.self.mean.mean.long, aes(x = factor(i), y = value), alpha = 1, size = 3) +  #stat = 'summary', fun.y = 'mean',
#  scale_y_continuous( breaks = c(seq(1, 3, by = 0.5), 4, 5, 10) ) +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Variance of species') +
  theme_bw()
p

## the scaling relation between mean and variance of individual species
slv1.XMeans.XVars.mean = inner_join(slv1.XMeans.mean.long, slv1.XVars.self.mean.long, by = c('i', 'variable'))

p  = ggplot(data = slv1.XMeans.XVars.mean, aes(x = value.x, y = value.y))
p  + geom_point() + 
  facet_wrap(~i, ncol=2, scales = 'free') + 
  scale_y_log10('log(species variance)') + scale_x_log10('log(species abundance)') +
  geom_smooth(method = 'lm', size = 1) +
  #  geom_text(data=eq, aes(x=9, y=1,  label=V1),parse = TRUE, inherit.aes=FALSE) +
  #  theme(strip.text.x = element_text(size = 14, angle = 0)) +
  theme_bw()

lm(log10(slv1.XMeans.mean[2,]) ~  log10(as.numeric(slv1.XVars.self.mean[2,1:50])))

## the summed variances and covariances of species
slv1.XVars.sum = aaply(slv1.XVars, .margins = c(1, 2), function(XVar) {
  c(varsum = sum(XVar), selfvarsum = sum(diag(XVar)), covsum = sum(XVar) - sum(diag(XVar)))
})
slv1.XVars.sum.mean = aaply(slv1.XVars.sum, .margins = c(1, 3), mean)
slv1.XVars.sum.mean = data.frame(slv1.XVars.sum.mean)
slv1.XVars.sum.mean$i = 1:dim(slv1.XVars)[1]
slv1.XVars.sum.mean.long = melt(slv1.XVars.sum.mean, id.vars = 'i')
p <- ggplot(slv1.XVars.sum.mean.long, aes(x = factor(i), y = value, group = variable, colour = variable)) +
  geom_point() + geom_line() +
  scale_color_discrete(breaks=c("varsum", "selfvarsum", "covsum"), name = '',
                       labels = c('summed variances and covariances', 'summed variances', 'summed covariances')) +
  theme(legend.title=element_blank()) +
  labs(x = 'netedness', y = 'summed variances and covariances') +
  theme_bw()
p


## temporal stability of individual species
slv1.si = slv1.XMeans.mean / sqrt(slv1.XVars.self.mean[, 1:50])
slv1.si = data.frame(slv1.si)
slv1.si$i = 1:dim(slv1.XMeans)[1]
slv1.si.long = melt(slv1.si, id.vars = 'i')
p <- ggplot(slv1.si.long, aes(x = factor(i), y = value, colour = factor(i))) +
  geom_point() + 
  geom_point(stat = 'summary', fun.y = 'mean', alpha = 1, size = 5) +  
  #scale_y_continuous( breaks = c(seq(1, 3, by = 0.5), 4, 5, 10) ) +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Temporal stability of individual species') +
  theme_bw()
p

## temporal stability of community
slv1.st = slv1.XMeans.mean.mean$slv1.XMeans.mean.mean / sqrt(slv1.XVars.sum.mean$varsum)
slv1.st = data.frame(slv1.st)
slv1.st$i = 1:dim(slv1.XMeans)[1]
p <- ggplot(slv1.st, aes(x = factor(i), y = slv1.st, group = 1, colour = factor(i))) +
  geom_point() + geom_line() +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Temporal stability of community') +
  theme_bw()
p  


###############################################################################
# check the deterministic LV1 model.
# species abundance in equilibrium <--> nested structure, competition and cooperation strengths
lv1.out = sim.lv1.2(graphs, beta0 = 0.001, gamma0 = 0.1)
lv1.Nstar = laply(lv1.out, function(lv1.one) {
  lv1.one$Nstar
})
lv1.Nstar = data.frame(lv1.Nstar)
lv1.Nstar$i = 1:124
lv1.Nstar = melt(lv1.Nstar, id.vars = 'i')
ggplot(lv1.Nstar, aes(x = factor(i), y = value, colour = factor(i))) +
  geom_point() + 
  geom_point(stat = 'summary', fun.y = 'mean', alpha = 1, size = 6) +
  guides(colour = FALSE ) +
  theme_bw()





XMeans = aaply(Xs, .margins = c(2,3), mean) # two dimensional matrix: (t+1)*m
matplot(XMeans, type = 'l', lwd = 0.7)

XVars = aaply(Xs, .margins = c(2), var) # three dimensional matrix: (t+1)*m*m

slv1.124.XVars.sum = aaply(slv1.124.XVars, .margins = c(1), function(XVar) {
  c(varsum = sum(XVar), selfvarsum = sum(diag(XVar)), covsum = sum(XVar) - sum(diag(XVar)))
})
matplot(slv1.124.XVars.sum, type = 'l', lwd = 0.7)

slv1.XVars.self = aaply(slv1.XVars, .margins = c(1,2), function(XVar) {
  c(diag(XVar))
})
matplot(XVars.self, type = 'l', lwd = 0.7)



XMeans = slv2_25_25_1.5_124.XMeans[,,502:1001,]
XMeans.mean = aaply(XMeans, .margins = c(1, 2), mean)
XMeans.sd = aaply(XMeans, .margins = c(1, 2), sd)

XVars = slv2_25_25_1.5_124.XVars
XVars.sum = aaply(XVars, .margins = c(1, 2, 3), function(XVar) {
  c(varsum = sum(XVar), selfvarsum = sum(diag(XVar)), covsum = sum(XVar) - sum(diag(XVar)))
})
matplot(XVars.sum[2,1,,], type = 'l', lwd = 0.7)

XVars.self = aaply(XVars, .margins = c(1, 2, 3), function(XVar) {
  c(diag(XVar))
})
matplot(XVars.self[2,1,,], type = 'l', lwd = 0.7)

XVars.cov = aaply(XVars, .margins = c(1, 2, 3), function(XVar) {
  XVar[lower.tri(XVar)]
})
matplot(XVars.cov[2,1,,1:100], type = 'l', lwd = 0.7)


###############################################################################
# the relation between the largest or smallest eigenvalues and the parms of matrix:
# 1) nested structure. 2) competition strength beta0. 3) cooperation strength gamma0
beta0s = seq(0, 1, by = 0.1)
gamma0s = seq(0, 0.2, by = 0.05)
A = ldply(beta0s, function(beta0) {
  ldply(gamma0s, function(gamma0) {
    ldply(graphs, function(graph) {
      graph.index = graph$count
      graph = graph$B
      edges = sum(graph > 0)
      numP = dim(graph)[1]
      numA = dim(graph)[2]
      D = diag(1, numP + numA)
      D[1:numP, 1:numP] = beta0
      D[(numP+1):(numP+numA), (numP+1):(numP+numA)] = beta0
      diag(D) = 1
      A = as.one.mode(graph)
      A[A > 0 ] = gamma0
      M = D - A  # competition interaction matrix
      c(beta0 = beta0, gamma0 = gamma0, graph = graph.index, levM = max(eigen(M)$values), 
        sevM = min(eigen(M)$values), levA = eigen(A)$values[1], levD = eigen(D)$values[1])  
    })
  })
})
ggplot(data = A, aes(x = graph, y = levM)) +
  geom_line() +
  facet_grid(beta0~gamma0)

A = sim.slv1.graph(graphs[[124]]$B, gamma0 = 0.148,steps = 1000, simnum = 500,sigma0 = 0.01)
B = sim.slv1.graph(graphs[[1]]$B, gamma0 = 0.148,steps = 1000, simnum = 500,sigma0 = 0.01)



slv1.0.1.XVars = aaply(slv1.out.0.1, .margins = c(2), var)
slv1.0.8.XVars = aaply(slv1.out.0.8, .margins = c(2), var)
mou.0.1.XVars = aaply(mou.out.0.1, .margins = c(2), var)
mou.0.8.XVars = aaply(mou.out.0.8, .margins = c(2), var)
mou.0.XVars = aaply(mou.out.0, .margins = c(2), var)

cummean2 <- function(A) {
  aaply(A, .margins = c(2), function(onecol) {
    cumsum(rev(onecol)) / seq_along(rev(onecol))
  })
}
cummean3 <- function(A) {
  aaply(A, .margins = c(1,3), function(onecol) {
    cumsum(rev(onecol)) / seq_along(rev(onecol))
  })
}
dim.reduct <- function(A) {
  aaply(A, .margins = c(1), function(oneslice) {
    c(diag(oneslice), oneslice[lower.tri(oneslice)])
  })
}

A = matrix(c(0.4465594,0.1578236,0.1578236, 0.3394577), ncol=2)
