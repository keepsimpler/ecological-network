library(plyr)
library(ggplot2)
library(dplyr)
library(reshape2)


## get temporal stability info from the mean vector and variance-covariance matrix
get.temporal.stability <- function(slv1.means.and.vars) {
  slv1.means = ldply(slv1.means.and.vars, function(one) {
    data.frame(means = one$means)
  })
  colnames(slv1.means) <- c('gindex', 'means')
  slv1.means.mean = slv1.means %.% group_by(gindex) %.% summarise(means.mean = mean(means))
  
  slv1.vars.self = ldply(slv1.means.and.vars, function(one) {
    data.frame(vars.self = diag(one$vars))
  })
  colnames(slv1.vars.self) <- c('gindex', 'vars.self')
  
  slv1.means.and.selfvars = slv1.means
  slv1.means.and.selfvars$vars.self = slv1.vars.self$vars.self

  slv1.means.and.selfvars.lm = dlply(slv1.means.and.selfvars, .(gindex), lm, formula = log10(means) ~ log10(vars.self))
  slv1.scaling = ldply(slv1.means.and.selfvars.lm, coef)
  colnames(slv1.scaling) <- c('gindex', 'intercept', 'power')
  
  slv1.means.and.selfvars$si = slv1.means.and.selfvars$means / sqrt(slv1.means.and.selfvars$vars.self)

  slv1.vars.sum = ldply(slv1.means.and.vars, function(one) {
    vars = one$vars
    c(varsum = sum(vars), selfvarsum = sum(diag(vars)), covsum = sum(vars) - sum(diag(vars)))
  })
  colnames(slv1.vars.sum)[1] <- 'gindex'
  slv1.vars.sum.long = melt(slv1.vars.sum, id.vars = c('gindex'))

  slv1.vars.sum$means.mean = slv1.means.mean$means.mean
  slv1.vars.sum$st = slv1.means.mean$means.mean / sqrt(slv1.vars.sum$varsum)
   
  list(means.selfvars.si = slv1.means.and.selfvars, meanofmeans.sumofvars.st = slv1.vars.sum )
}


slv1.means = ldply(slv1.means.and.vars, function(one) {
  data.frame(means = one$means)
})
#slv1.means$i = rep(1:18, each = 20)
colnames(slv1.means) <- c('gindex', 'means')
slv1.means.mean = slv1.means %.% group_by(gindex) %.% summarise(means.mean = mean(means))

#slv1.means = slv1.means %.% filter(i < 35)
p <- ggplot(slv1.means, aes(x = factor(gindex), y = means, colour = factor(gindex))) +
  geom_point() + 
  geom_point(data = slv1.means.mean, aes(x = factor(gindex), y = means.mean), alpha = 1, size = 5) +  
#  scale_y_continuous( breaks = c(seq(1, 3, by = 0.5), 4, 5, 10) ) +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Biomass of species') +
  theme_bw()
p


slv1.vars.self = ldply(slv1.means.and.vars, function(one) {
  data.frame(vars.self = diag(one$vars))
})
#slv1.vars.self$i = rep(1:36, each = 20)
colnames(slv1.vars.self) <- c('gindex', 'vars.self')

slv1.means.and.selfvars = slv1.means
slv1.means.and.selfvars$vars.self = slv1.vars.self$vars.self
p  = ggplot(data = slv1.means.and.selfvars, aes(x = means, y = vars.self))
p  + geom_point() + 
  facet_wrap(~ gindex, ncol = 6, scales = 'free') + 
  scale_y_log10('log(species variance)') + scale_x_log10('log(species abundance)') +
  geom_smooth(method = 'lm', size = 1) +
  #  geom_text(data=eq, aes(x=9, y=1,  label=V1),parse = TRUE, inherit.aes=FALSE) +
  #  theme(strip.text.x = element_text(size = 14, angle = 0)) +
  theme_bw()

slv1.means.and.selfvars.lm = dlply(slv1.means.and.selfvars, .(gindex), lm, formula = log10(means) ~ log10(vars.self))
slv1.scaling = ldply(slv1.means.and.selfvars.lm, coef)
colnames(slv1.scaling) <- c('gindex', 'intercept', 'power')

slv1.means.and.selfvars$si = slv1.means.and.selfvars$means / sqrt(slv1.means.and.selfvars$vars.self)
#slv1.means = slv1.means %.% filter(i < 35)
p <- ggplot(slv1.means.and.selfvars, aes(x = factor(gindex), y = si, colour = factor(gindex))) +
  geom_point() + 
#  geom_point(data = slv1.means.mean, aes(x = factor(gindex), y = means.mean), alpha = 1, size = 5) +  
  #  scale_y_continuous( breaks = c(seq(1, 3, by = 0.5), 4, 5, 10) ) +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Temporal stability of individual species') +
  theme_bw()
p


slv1.vars.sum = ldply(slv1.means.and.vars, function(one) {
  vars = one$vars
  c(varsum = sum(vars), selfvarsum = sum(diag(vars)), covsum = sum(vars) - sum(diag(vars)))
})
colnames(slv1.vars.sum)[1] <- 'gindex'
slv1.vars.sum.long = melt(slv1.vars.sum, id.vars = c('gindex'))
p <- ggplot(slv1.vars.sum.long, aes(x = factor(gindex), y = value, group = variable, colour = variable)) +
  geom_point() + geom_line() +
  scale_color_discrete(breaks=c("varsum", "selfvarsum", "covsum"), name = '',
                       labels = c('summed variances and covariances', 'summed variances', 'summed covariances')) +
  theme(legend.title=element_blank()) +
  labs(x = 'netedness', y = 'summed variances and covariances') +
  theme_bw()
p

slv1.st = data.frame(gindex = slv1.means.mean$gindex)
slv1.st$st = slv1.means.mean$means.mean / sqrt(slv1.vars.sum$varsum)
p <- ggplot(slv1.st, aes(x = factor(gindex), y = st, group = 1, colour = factor(gindex))) +
  geom_point() + geom_line() +
  guides(colour = FALSE ) +
  labs(x = 'netedness', y = 'Temporal stability of community') +
  theme_bw()
p  
