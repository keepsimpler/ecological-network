#!/usr/bin/Rscript
library(plyr)

library(doMC)  # 
registerDoMC()  # register Multi Cores
getDoParWorkers()  # get available Cores

load(file = 'slv1.20.RData')
slv1.means.and.vars = alply(A,.parallel = TRUE, .margins = c(1), function(one) {
  print(1)
  one.means = aaply(one, .margins = c(2, 3), mean)
  one.means = one.means[5002:10001, ]
  one.means.mean = aaply(one.means, .margins = c(2), mean)
  
  one.vars = aaply(one, .margins = c(2), var)
  one.vars = one.vars[5002:10001, , ]
  one.vars.mean = aaply(one.vars, .margins = c(2, 3), mean)
  
  list(means = one.means.mean, vars = one.vars.mean)
})

slv1.means.and.vars[[length(slv1.means.and.vars) + 1]] = list(gamma)

save(slv1.means.and.vars, file = 'slv1.means.and.vars.RData')