library(plyr)
library(reshape2)
library(ggplot2)

### Estimate the largest eigenvalues of empirical mutualistic networks

df.wigner.quan =   # using semicircle plus twins method
  ldply(datasets.connected.pre, function(A) {
    A = get(A)
    est.lev.wigner(A, is.binary = FALSE)
  })
df.wigner.quan['dataset'] = datasets.connected.pre
df.wigner.quan$error.wigner = abs((df.wigner.quan$lev.est - df.wigner.quan$lev.true) / df.wigner.quan$lev.true)
df.wigner.quan$error.chung = abs((df.wigner.quan$chung - df.wigner.quan$lev.true) / df.wigner.quan$lev.true)

df.semisuperellipse.quan =   # using semisuperellipse plus twins method
  ldply(datasets.connected.pre, function(A) {
    A = get(A)
    est.lev.semisuperellipse(A, is.binary = FALSE)
  })
df.semisuperellipse.quan['dataset'] = datasets.connected.pre
df.semisuperellipse.quan$error.semisuperellipse = 
  abs((df.semisuperellipse.quan$lev.est - df.semisuperellipse.quan$lev.true) / df.semisuperellipse.quan$lev.true)

df.errors = df.wigner.quan[c('dataset', 'error.chung', 'error.wigner')]
df.errors['error.semisuperellipse'] = df.semisuperellipse.quan$error.semisuperellipse
df.errors['type'] = 'Quantitative'

df.wigner.bin =   # using semicircle plus twins method
  ldply(datasets.connected.pre, function(A) {
    A = get(A)
    est.lev.wigner(A, is.binary = TRUE)
  })
df.wigner.bin['dataset'] = datasets.connected.pre
df.wigner.bin$error.wigner = abs((df.wigner.bin$lev.est - df.wigner.bin$lev.true) / df.wigner.bin$lev.true)
df.wigner.bin$error.chung = abs((df.wigner.bin$chung - df.wigner.bin$lev.true) / df.wigner.bin$lev.true)

df.semisuperellipse.bin =   # using semisuperellipse plus twins method
  ldply(datasets.connected.pre, function(A) {
    A = get(A)
    est.lev.semisuperellipse(A, is.binary = TRUE)
  })
df.semisuperellipse.bin['dataset'] = datasets.connected.pre
df.semisuperellipse.bin$error.semisuperellipse = 
  abs((df.semisuperellipse.bin$lev.est - df.semisuperellipse.bin$lev.true) / df.semisuperellipse.bin$lev.true)

tmp = df.wigner.bin[c('dataset', 'error.chung', 'error.wigner')]
tmp['error.semisuperellipse'] = df.semisuperellipse.bin$error.semisuperellipse
tmp['type'] = 'Binary'

df.errors = rbind(df.errors, tmp)
df.errors.long = melt(df.errors, id.vars = c('dataset', 'type'), measure.vars = c('error.wigner', 'error.chung', 'error.semisuperellipse'),
     variable.name = 'error.type', value.name = 'error.value')


p <- ggplot(data=df.errors.long, aes(x=error.type,y=error.value))
p + geom_boxplot() +
  facet_wrap(~type, scale='free') +
  scale_x_discrete('Estimation Methods', labels=c('Semicircle plus \ntwins', 'Chung',
                                                  'Semi-superellipse \nplus twins')) +
  scale_y_continuous('|Relative Errors|') +
  #  theme(axis.text.x = element_text(angle = 90, hjust = 0)) + 
  theme_bw()

