################################################################################
# Analyze and display the results

library(plyr)
library(reshape2)
library(ggplot2)

## Input: A four dimensional array, 
## First dimension: the index of graphs with increasing nestedness
## Second dimension: the number of simulation of Stochastic LV2 model
## Third dimension: the time steps of simulation
## Fourth dimension: the number of species
A = slv1.XMeans[, 2:1001, ]
B = slv1.XVars.self[, 2:1001, ]
C = A / sqrt(B)
display.1 <- function(A) {
  A = A[, 501:1000, ]  # subset the stationary data
  B = aaply(A, .margins = c(1, 3), mean)
  B = data.frame(B)
  B$i = graphs.index
  C = melt(B, id.vars = 'i')
  p <- ggplot(C, aes(x = factor(i), y = value, colour = factor(i))) +
    geom_point() + 
    geom_point(stat = 'summary', fun.y = 'mean', alpha = 1, size = 6) +
    guides(colour = FALSE ) +
    theme_bw()
  p
}

A = slv1.XVars.sum[, 2:1001, ]
A = A[, 501:1000, ]
B = aaply(A, .margins = c(1, 3), mean)
matplot(B, type = 'l', lwd = 0.7)


## check the deterministic LV1
lv1.out = sim.lv1(graphs)
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

slv1.XVars.sum = aaply(slv1.XVars, .margins = c(1,2), function(XVar) {
  c(varsum = sum(XVar), selfvarsum = sum(diag(XVar)), covsum = sum(XVar) - sum(diag(XVar)))
})
matplot(XVars.sum, type = 'l', lwd = 0.7)

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
