###############################################################################
## Estimate the largest eigenvalues of empirical mutualistic networks

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
df.errors.long = melt(df.errors, id.vars = c('dataset', 'type'), measure.vars = c('error.semisuperellipse', 'error.wigner', 'error.chung'),
     variable.name = 'error.type', value.name = 'error.value')

# the accuray and precision of three estimate method for both binary and quantitative bipartite networks
df.errors.accuracy = ddply(df.errors.long, error.type~type, summarise, 
                           accuracy = mean(error.value), precision = sd(error.value))

# Figure 3 in main text
p <- ggplot(data=df.errors.long, aes(x=error.type,y=error.value))
p + geom_boxplot() +
  facet_wrap(~type, scale='free') +
  scale_x_discrete('Estimation Methods', labels=c('Semi-superellipse \nplus twins', 'Semicircle plus \ntwins', 'Chung')) +
  scale_y_continuous('|Relative Errors|') +
  geom_text(data = df.errors.accuracy, aes(x = error.type, y = 0.7,
                                           label = paste(sprintf('%.4f',accuracy), '+-', sprintf('%.4f',precision))), parse = TRUE) +
  theme_bw()



################################## Null Models #############################################
#library(doMC)
#registerDoMC()
#getDoParWorkers()

# measures of empirical networks
df.empirical = ldply(datasets.connected.pre, .parallel = TRUE, function(A) {
  A = get(A)
  get.measures(A)
})
df.empirical$dataset = datasets.connected.pre


tried = 1:1000  # number of randomizations

########### measures of Null Model 1 #############
df.nm1 = ldply(datasets.connected.pre, .parallel = TRUE, function(A) {
  A = get(A)
  ldply(tried, .parallel = TRUE, function(i) {
    print(paste(i))
    B = A
    B[B > 0] <- sample(A[A > 0])  # NM1 : shuffling weights among the existing links
    get.measures(B)
  })
})
df.nm1$tried = rep(tried, times = length(datasets.connected.pre))
df.nm1$dataset = rep(datasets.connected.pre, each = length(tried))

df.nm1 = merge(df.nm1, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm1.pvalue = ddply(df.nm1, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / length(tried), 
                      zvalue.lev = 1- pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)),
                      pvalue.Hq = sum(Hq.x > Hq.y) / length(tried), 
                      zvalue.Hq = 1- pnorm((mean(Hq.y) - mean(Hq.x)) / sd(Hq.x))
                      )
sum(df.nm1.pvalue$pvalue.lev > 0.95 | df.nm1.pvalue$pvalue.lev < 0.05)  # pvalues not in [0.05, 0.95]
sum(df.nm1.pvalue$pvalue.Hq > 0.95 | df.nm1.pvalue$pvalue.Hq < 0.05)  # pvalues not in [0.05, 0.95]

# get the correlation between Hq and lev
df.nm1.corr = ddply(df.nm1, .(dataset), summarise, pearson = cor(sqrt(sqrt(Hq.x)), lev.weight.x),  # peason coefficient between Hq and lev
      R2 = summary(lm(sqrt(sqrt(Hq.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between Hq and lev
      label.position.y = min(lev.weight.x) + 0.1 * (max(lev.weight.x) - min(lev.weight.x)),  # plot label position
      label.position.x =  (min(sqrt(sqrt(Hq.x))) + 0.7 * (max(sqrt(sqrt(Hq.x))) - min(sqrt(sqrt(Hq.x))))))  

# plot the correlation between Hq and lev
ggplot(data = df.nm1, mapping = aes(x = sqrt(sqrt(Hq.x)), y = lev.weight.x)) + 
  geom_point(size=1, alpha=.7, color='sky blue', fill='sky blue') + 
  geom_point(data = df.nm1, mapping = aes(x = sqrt(sqrt(Hq.y)), y = lev.weight.y),
             size=2.5, alpha=1, color='red') + 
  facet_wrap(~ dataset, ncol = 5, scales = 'free') +
  geom_smooth(method = 'lm', size = 0.5, linetype='dashed', color = 'blue') +
  theme_bw() +
  scale_x_continuous(expression(sqrt(H[q],4))) +
  scale_y_continuous(expression(lambda[1]^nm1)) +
  geom_text(data = df.nm1.corr, aes(x = label.position.x, y = label.position.y,
                                                      label = paste(expression(R^2), '==', sprintf('%.2f',R2))), parse = T) +
  theme(axis.title = element_text(size = rel(1.4), face = 'bold'), 
        strip.text.x = element_text(size = rel(1.35)))

# plot the pvalues of lev
ggplot(data = df.nm1, mapping = aes(x = lev.weight.x / lev.weight.y, y = ..scaled.., ymin=0)) + 
  geom_density(size = 1, alpha = .3, color = 'sky blue', fill = 'sky blue') + 
  geom_vline(aes(xintercept = 1), linetype = 'dashed', size = 1, color = 'red')  +
  facet_wrap(~ dataset, ncol=5, scales = 'free_y') +
  theme_bw() + coord_equal() +
  scale_x_continuous(expression(lambda[1]^NM1/lambda[1])) +
  scale_y_continuous('Density') +
  theme(axis.title=element_text(size=14, face='bold')) +
  geom_text(data=df.nm1.pvalue, aes(1.15, 0.9, label=paste('p =', pvalue.lev), color='red')) +
  theme(legend.position = "none")

p = ggplot(data = df.nm1.pvalue, aes(x = dataset, y = pvalue.lev))
p + geom_point() + theme_bw() +
  geom_hline(yintercept=0.05, linetype='dashed', size=1, color='red')  +
  geom_hline(yintercept=0.95, linetype='dashed', size=1, color='red')  +
  scale_x_discrete('datasets', labels=c(1:40)) +
  scale_y_continuous('p-values', breaks=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00))


########### measures of Null Model 2 #############
sigmas = seq(0.1,1.0,0.1)  # the parameter of log-normal distribution
df.nm2 =   
  ldply(datasets.connected.pre, .parallel = FALSE, function(A) {  # for every dataset
    A = get(A)
    totallinks <- length(A[A > 0])  # get number of links
    totalweights <- sum(A)  # get sum of link weights
    ldply(sigmas, .parallel = TRUE, function(sigma) {  # for every \sigma
      print(paste(sigma))
      ldply(tried, .parallel = TRUE, function(i) {  # try [tried] times
        B = A
        #  NM2: randomizing link weights according to log-normal distribution, keep total weights fixed
        B[B > 0] = rlnorm.fixsum(n = totallinks, total = totalweights, sigmalog = sigma)
        get.measures(B)
        })
    })
  })
df.nm2$tried = rep(tried, times = length(sigmas) * length(datasets.connected.pre))
df.nm2$sigmas = rep(rep(sigmas, each = length(tried)), times = length(datasets.connected.pre))
df.nm2$dataset = rep(datasets.connected.pre, each = length(sigmas) * length(tried))

tmp = merge(df.nm2, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm2.pvalue = ddply(tmp, .(dataset, sigmas), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / length(tried), 
                      zvalue.lev = 1- pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)),
                      pvalue.W2 = sum(w2.x > w2.y) / length(tried), 
                      zvalue.W2 = 1- pnorm((mean(w2.y) - mean(w2.x)) / sd(w2.x)),
                      pvalue.Hq = sum(Hq.x > Hq.y) / length(tried), 
                      zvalue.Hq = 1- pnorm((mean(Hq.y) - mean(Hq.x)) / sd(Hq.x)))

p <- ggplot(data = df.nm2.pvalue, mapping = aes(x = sigmas))
p +  geom_line(aes(y=pvalue.lev, colour="lambda1")) + 
  geom_line(aes(y=pvalue.Hq, colour="Hq")) + 
  geom_line(aes(y=pvalue.W2, colour="W2")) + 
  #  geom_line(aes(y=pvalue.nest,  colour="WNODF"),linetype='dashed') + 
  #  geom_line(aes(y=pvalue.wine, colour="WINE"), linetype='dotted') + 
  scale_colour_manual("",  breaks = c("lambda1", "Hq", "W2"), #, "WNODF", "WINE" 
                      values = c("red", "green","blue")) +  #, "grey","grey"
  scale_x_continuous(expression(sigma)) +
  scale_y_continuous('p-values') +
  facet_wrap(~dataset, ncol=5, scales = 'free') + 
  theme_bw() +
  theme(legend.key = element_rect(colour = "white"), 
        axis.title = element_text(size = rel(1.3), face = 'bold'),
        strip.text.x = element_text(size = rel(1.1)))


########### measures of Null Model 3 #############
df.nm3 = ldply(datasets.connected.pre, .parallel = TRUE, function(A) {
  A = get(A)
  ldply(tried, .parallel = TRUE, function(i) {
    print(paste(i))
    B = swaplinks(A, HowManyToTry = 2000)  # NM3: swaping links that keep the degree distribution constant
    get.measures(B)
  })
})
df.nm3$tried = rep(tried, times = length(datasets.connected.pre))
df.nm3$dataset = rep(datasets.connected.pre, each = length(tried))

df.nm3.random = df.nm3  # backup 
df.nm3 = merge(df.nm3, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm3.pvalue = ddply(df.nm3, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / length(tried), 
                      zvalue.lev = 1 - pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)),
                      pvalue.Hq = sum(Hq.x > Hq.y) / length(tried), 
                      zvalue.Hq = 1- pnorm((mean(Hq.y) - mean(Hq.x)) / sd(Hq.x)),
                      pvalue.NODF = sum(NODF.x > NODF.y) / length(tried), 
                      zvalue.NODF = 1- pnorm((mean(NODF.y) - mean(NODF.x)) / sd(NODF.x)),
                      pvalue.CMNB = sum(cmnb.x > cmnb.y) / length(tried), 
                      zvalue.CMNB = 1- pnorm((mean(cmnb.y) - mean(cmnb.x)) / sd(cmnb.x)),
                      pvalue.WNODF = sum(WNODF.x > WNODF.y) / length(tried), 
                      zvalue.WNODF = 1- pnorm((mean(WNODF.y) - mean(WNODF.x)) / sd(WNODF.x)),
                      pvalue.WINE = sum(wine.x > wine.y) / length(tried), 
                      zvalue.WINE = 1- pnorm((mean(wine.y) - mean(wine.x)) / sd(wine.x))
)
sum(df.nm3.pvalue$pvalue.CMNB > 0.95 | df.nm3.pvalue$pvalue.CMNB < 0.05)  # pvalues not in [0.05, 0.95]

# get the correlation between Hq and lev
colnames(df.nm3)[16] = 'WNODF.x'; colnames(df.nm3)[33] = 'WNODF.y'  # rename the column name of 'weighted NODF' to 'WNODF'
df.nm3.corr = ddply(df.nm3, .(dataset), summarise, pearson.Hq = cor(sqrt(sqrt(Hq.x)), lev.weight.x),  # peason coefficient between Hq and lev
                    R2.Hq = summary(lm(sqrt(sqrt(Hq.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between Hq and lev
                    pearson.NODF = cor(sqrt(sqrt(NODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.NODF = summary(lm(sqrt(sqrt(NODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.CMNB = cor(sqrt(sqrt(cmnb.x)), lev.weight.x),  # peason coefficient between cmnb and lev
                    R2.CMNB = summary(lm(sqrt(sqrt(cmnb.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between cmnb and lev
                    pearson.WNODF = cor(sqrt(sqrt(WNODF.x)), lev.weight.x),  # peason coefficient between WNODF and lev
                    R2.WNODF = summary(lm(sqrt(sqrt(WNODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between WNODF and lev
                    pearson.WINE = cor(sqrt(sqrt(wine.x)), lev.weight.x),  # peason coefficient between WINE and lev
                    R2.WINE = summary(lm(sqrt(sqrt(wine.x)) ~ lev.weight.x, data=faithful))$r.squared)  # R^2 coefficient between WINE and lev
                    
df.nm3.corr.long = melt(df.nm3.corr, id.vars = c('dataset'), measure.vars = c('R2.Hq', 'R2.NODF', 'R2.CMNB', 'R2.WNODF', 'R2.WINE'),
                      variable.name = 'corr.type', value.name = 'corr.value')
df.nm3.corr.long$nullmodel = 'Null model 3'



########### measures of Null Model 4 #############
df.nm4 = ldply(datasets.connected.pre, .parallel = TRUE, function(A) {
  A = get(A)
  edges = length(A[A > 0])  # get number of links
  nodes = dim(A)  # number of nodes of two groups
  ldply(tried, .parallel = TRUE, function(i) {
    print(paste(i))
    # NM4: generate  a ER random bipartite graph with same number of nodes and edges with empirical networks
    g = bipartite.random.game(nodes[1], nodes[2], type = 'gnm', m = edges)
    B = get.incidence(g)  # get the incidence matrix of the random bipartite graph
    B[B > 0] = A[A > 0]  # keep the link weights of the empirical network
    get.measures(B)  # get measures of the randomized network
  })
})
df.nm4$tried = rep(tried, times = length(datasets.connected.pre))
df.nm4$dataset = rep(datasets.connected.pre, each = length(tried))

df.nm4.random = df.nm4  # backup 
df.nm4 = merge(df.nm4, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm4.pvalue = ddply(df.nm4, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / length(tried), 
                      zvalue.lev = 1 - pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)))
sum(df.nm4.pvalue$pvalue.lev > 0.95 | df.nm4.pvalue$pvalue.lev < 0.05)  # pvalues not in [0.05, 0.95]

# get the correlation between Hq and lev
colnames(df.nm4)[16] = 'WNODF.x'; colnames(df.nm4)[33] = 'WNODF.y'  # rename the column name of 'weighted NODF' to 'WNODF'
df.nm4.corr = ddply(df.nm4, .(dataset), summarise, pearson.Hq = cor(sqrt(sqrt(Hq.x)), lev.weight.x),  # peason coefficient between Hq and lev
                    R2.Hq = summary(lm(sqrt(sqrt(Hq.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between Hq and lev
                    pearson.NODF = cor(sqrt(sqrt(NODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.NODF = summary(lm(sqrt(sqrt(NODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.CMNB = cor(sqrt(sqrt(cmnb.x)), lev.weight.x),  # peason coefficient between cmnb and lev
                    R2.CMNB = summary(lm(sqrt(sqrt(cmnb.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between cmnb and lev
                    pearson.WNODF = cor(sqrt(sqrt(WNODF.x)), lev.weight.x),  # peason coefficient between WNODF and lev
                    R2.WNODF = summary(lm(sqrt(sqrt(WNODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between WNODF and lev
                    pearson.WINE = cor(sqrt(sqrt(wine.x)), lev.weight.x),  # peason coefficient between WINE and lev
                    R2.WINE = summary(lm(sqrt(sqrt(wine.x)) ~ lev.weight.x, data=faithful))$r.squared)  # R^2 coefficient between WINE and lev

df.nm4.corr.long = melt(df.nm4.corr, id.vars = c('dataset'), measure.vars = c('R2.Hq', 'R2.NODF', 'R2.CMNB', 'R2.WNODF', 'R2.WINE'),
                        variable.name = 'corr.type', value.name = 'corr.value')
df.nm4.corr.long$nullmodel = 'Null model 4'

df.nm34.corr.long = rbind(df.nm3.corr.long, df.nm4.corr.long)

p <- ggplot(data=df.nm34.corr.long, aes(x = corr.type, y = corr.value))
p + geom_boxplot() +
    facet_wrap(~nullmodel, scale='free') +
  scale_x_discrete('Network measure', labels = c('Hq', 'NODF', 'CMNB', 'WNODF', 'WINE')) +
  scale_y_continuous(expression(paste('Coefficient of determination ', R^2))) +
  theme_bw()




##########################################################################################
## Plot Semicircle and Semi-superellipse Plus Twins hyperthesis for Bipartite networks
##########################################################################################

####################################################################################################
## the largest eigenvalue grows faster than the second eigenvalue for ER random graphs. 
ns = seq(100, 1000, by = 100)
gap.er = ldply(ns, .parallel = TRUE, function(n) {
  edgess = seq(n * 3, n * 15, by = n)
  ldply(edgess, .parallel = TRUE, function(edges) {
    print(paste(n, edges))
    G = erdos.renyi.game(n, edges, type="gnm")
    A = igraph::get.adjacency(G)
    lev = eigen(A)$values[1]
    sev = eigen(A)$values[2]
    mev = eigen(A)$values[n]
    r = (lev - sev) / (sev - mev)
    c(n, edges, r)
  })
})
colnames(gap.er) <- c('n', 'edges', 'r')
ggplot(data = gap.er, aes(x = edges, y = r)) +
  geom_line() +
  facet_wrap(~n, scales = 'free')

####################################################################################################
## the largest eigenvalue grows faster than the second eigenvalue for Scale-Free random graphs
ns = seq(100, 1000, by = 100)
powers = seq(3.1, 4.5, by = 0.1)
gap.sf = ldply(ns, function(n) {
  edgess = seq(n * 3, n * 15, by = n)
  ldply(edgess, .parallel = TRUE, function(edges) {
    print(paste(n, edges))    
    ldply(powers, .parallel = TRUE, function(power) {
      G = static.power.law.game(n, edges, power)
      A = igraph::get.adjacency(G)
      lev = eigen(A)$values[1]
      sev = eigen(A)$values[2]
      mev = eigen(A)$values[n]
      r = (lev - sev) / (sev - mev)
      c(n, edges, power, r)
    })
  })
})
colnames(gap.sf) <- c('n', 'edges', 'power', 'r')
gap.sf2 = ddply(gap.sf, .(n, power), summarise, r = mean(r))
ggplot(data = gap.sf2, aes(x = power, y = r)) +
  geom_line() +
  facet_wrap(~n, scales = 'free')

####################################################################################################
## the largest eigenvalue grows faster than the second eigenvalue for bipartite ER random graphs.
ns = seq(50, 500, by = 50)  # a sequence of number of nodes
ks = seq(3, 15, by = 1)  # a sequence of average node degree
tried = 1:100  # number of randomizations for each n and k
gap.bipartite.er = ldply(ns, .parallel = FALSE, function(n) {
  ldply(ks, .parallel = TRUE, function(k) {
    print(paste(n, k))
    ldply(tried, .parallel = TRUE, function(i) {
      G = graph.connected(s = c(n, n), k = k, gtype = 'bipartite')  # generate a random connected bipartite graph
      A = igraph::get.adjacency(G)
      A = as.matrix(A)
      lev = eigen(A)$values[1]
      sev = eigen(A)$values[2]
      r = lev / sev
      c(n, k, i, r, lev, sev)
    })
  })
})
# gap.bipartite.er$c = gap.bipartite.er$k / gap.bipartite.er$n
colnames(gap.bipartite.er) <- c('n', 'k', 'i', 'r', 'lev', 'sev')
gap.bipartite.er2 = ddply(gap.bipartite.er, .(n, k), summarise, r = mean(r) )
library(akima)
im = with(gap.bipartite.er2, interp(n, k, r))  
with(im,
     filled.contour(x, y, z, col = rainbow(200), nlevels = 200,
                    plot.title = title(main = 'The spectral gap of bipartite ER random graphs', 
                                       xlab = 'Number of species (s)', ylab = 'Average species degree (<k>)'), 
                    key.title = title(main = expression(paste(lambda[1]/lambda[2], sep = ''))),
                    key.axes = axis(4, seq(1.5, 5, by = 0.5)),
                    plot.axes = {axis(1); axis(2) ; contour(x,y,z, add = T, size = 15); contour(x,y,z, nlevels = 1, level = 0, add = T, lwd = 1.5)}
     )
)


###############################################################################
## Example of the semi-superellipse plus twins law of bipartite ER random graph
###############################################################################
semisuperellipse.plus.twins.of.bipartite <- function(k, n, title) {
  k = k  # 2 5 10
  G = graph.connected(s = c(100, 100), k = k, gtype = 'bipartite')  # generate a random connected bipartite graph
  A = igraph::get.adjacency(G)  # get 
  A = as.matrix(A)
  s = dim(A)[1]
  lev = eigen(A)$values[1]
  sev = eigen(A)$values[2]
  eigenvalues = eigen(A)$values
  H <- hist(eigenvalues, breaks = 75, plot = FALSE) 
  Hbreaks <- H$breaks
  Hmids <- H$mids
  Hcounts <- H$counts 
  temp <- data.frame(
    k = k,
    x = Hmids,
    y = Hcounts
  )
  y = FindSECounts(Hbreaks, Hbreaks[-1], r = sev, n = n) * (s - 1)  # n = 0.8  1.5 2.0
  temp2 <- data.frame(
    k = k,
    x = Hmids,
    y = y
  )
  
  p <- ggplot(data = temp, aes(x = x, ymax = y, ymin = 0)) +  # p1 p2 p3
    geom_linerange(size = 0.75) +
    geom_line(data = temp2, aes(x = x, y = y), colour = "red") +
    scale_x_continuous(expression(paste(lambda))) +  # (k = 2 5 10)
    scale_y_continuous('Density') +
    ggtitle(title) +  # 'a) (k = 2)'
    theme_bw()
 p 
}
p1 <- semisuperellipse.plus.twins.of.bipartite(k = 2, n = 0.8, title = '<k> = 2')
p2 <- semisuperellipse.plus.twins.of.bipartite(k = 5, n = 1.5, title = '<k> = 5')
p3 <- semisuperellipse.plus.twins.of.bipartite(k = 10, n = 2.0, title = '<k> = 10')

multiplot(p1, p2, p3, cols=3)



###############################################################################
## Example of null models
###############################################################################
library(corrplot)  # 
library(igraph)
col3 <- colorRampPalette(c("red", "white", "blue"))
A = olesen2002flores.C.Pre  # dataset
# A = sortweb(A>0)
numP = dim(A)[1]  # plants
numA = dim(A)[2]  # animals
B = matrix(sample(x = c(-1,1), replace = TRUE, size = numP * numA), ncol = numA, nrow = numP)
B = A * B  # for more colors, randomly change values of some elements to their negitive values
## The empirical network
corrplot(B, method = "color", cl.pos = "n", tl.pos = "n", addgrid.col = 'black', 
#         title = 'Empirical ecological network',
         is.corr = FALSE, p.mat = abs(B),  insig = "p-value", 
         col = col3(200))  #, cl.lim = c(min(B),max(B))

## Null Model 1
NM1 = B
NM1[NM1 != 0] <- sample(B[B != 0]) 
corrplot(NM1, method = "color", cl.pos = "n", tl.pos = "n", addgrid.col = 'black', 
         is.corr = FALSE, p.mat = abs(NM1),  insig = "p-value", 
         col = col3(200))

## Null Model 2
NM2 = B
totallinks <- length(B[B != 0])  # get number of links
totalweights <- sum(abs(A))  # get sum of link weights
#  randomizing link weights according to log-normal distribution, keep total weights fixed
NM2[NM2 != 0] = rlnorm.fixsum(n = totallinks, total = totalweights, sigmalog = 0.1)
NM2 = NM2 *  # for more colors, 
  matrix(sample(x = c(-1,1), replace = TRUE, size = numP * numA), ncol = numA, nrow = numP)
corrplot(NM2, method = "color", cl.pos = "n", tl.pos = "n", addgrid.col = 'black', 
         is.corr = FALSE, p.mat = abs(NM2),  insig = "p-value", 
         col = col3(200))

## Null Model 3
NM3 = swaplinks(B, HowManyToTry = 5000)  # swaping links that keep the degree distribution constant
corrplot(NM3, method = "color", cl.pos = "n", tl.pos = "n", addgrid.col = 'black', 
         is.corr = FALSE, p.mat = abs(NM3),  insig = "p-value", 
         col = col3(200))

## Null Model 4
g = bipartite.random.game(numP, numA, type = 'gnm', m = totallinks)
NM4 = get.incidence(g)  # get the incidence matrix of the random bipartite graph
NM4[NM4 != 0] = B[B != 0]  # keep the link weights of the empirical network
corrplot(NM4, method = "color", cl.pos = "n", tl.pos = "n", addgrid.col = 'black', 
         is.corr = FALSE, p.mat = abs(NM4),  insig = "p-value", 
         col = col3(200))

