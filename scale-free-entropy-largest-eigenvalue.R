require('igraph')
require('ggplot2')
n = 200  # node number
d = 5 # node degree
edges = n * d #
tried = 500

maxeig.entropy.sf.nm = data.frame(alpha=numeric(0), maxeig.sf=numeric(0), entropy.sf=numeric(0), squaresum.degrees.sf=numeric(0))
for (i in seq(from=2, to=6, by=0.01)) {
  #gnm = erdos.renyi.game(n, edges, type="gnm")
  repeat {
    gsf = static.power.law.game(n, edges, i)
    if (is.connected(gsf)) break
  }
  asf = get.adjacency(gsf)
  maxeig.sf = max(Re(eigen(asf)$values))
  degrees.sf = degree(gsf)
  degrees.sf = degrees.sf / sum(degrees.sf)
  squaresum.degrees.sf = sum(degrees.sf^2)
  log.degrees.sf = degrees.sf * log(degrees.sf)
  entropy.sf = - sum(log.degrees.sf)

  newrow = c(i, maxeig.sf, entropy.sf, squaresum.degrees.sf)
  maxeig.entropy.sf.nm = rbind(maxeig.entropy.sf.nm, newrow)
}
colnames(maxeig.entropy.sf.nm) = c('alpha', 'maxeig.sf', 'entropy.sf', 'squaresum.degrees.sf')
ggplot(data=maxeig.entropy.sf.nm, aes(x=entropy.sf, y=maxeig.sf)) +
  geom_line() +
  geom_smooth(method = 'lm')

lm(maxeig.sf ~ entropy.sf, data = maxeig.entropy.sf.nm)

gnm = static.power.law.game(n, edges, 3)
#gnm = erdos.renyi.game(n, edges, type="gnm",directed = F)
anm = get.adjacency(gnm)
anm = as.matrix(anm)
estimate.dominant.eigenvalue(anm)

maxeig.sf = max(Re(eigen(anm)$values))
banm = rbind(cbind(matrix(rep(0, n * n), ncol = n), as.matrix(anm)),
             cbind(t(as.matrix(anm)), matrix(rep(0, n * n), ncol = n)))
hist(Re(eigen(banm)$values))
# s = 500
# mean = -1
# std = 2
# M = matrix(rnorm(s*s, mean, std),ncol=s)
# eigs = eigen(M)$values
# max(Re(eigs))
# P = degree.animals %o% degree.plants
