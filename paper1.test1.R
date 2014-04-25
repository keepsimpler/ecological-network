library(bipartite)
library(igraph)
library(plyr)
library(reshape2)
library(ggplot2)

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
library(doMC)

registerDoMC()
getDoParWorkers()

# measures of empirical networks
df.empirical = ldply(datasets.connected.pre, function(A) {
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
df.nm1.pvalue = ddply(df.nm1, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / 500, 
                      zvalue.lev = 1- pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)))
sum(df.nm1.pvalue$pvalue.lev > 0.95 | df.nm1.pvalue$pvalue.lev < 0.05)  # pvalues not in [0.05, 0.95]

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
      ldply(tried, .parallel = TRUE, function(i) {  # try 500 times
        B = A
        #  randomizing link weights according to log-normal distribution, keep total weights fixed
        B[B > 0] = rlnorm.fixsum(n = totallinks, total = totalweights, sigmalog = sigma)
        get.measures(B)
        })
    })
  })
df.nm2$tried = rep(tried, times = length(sigmas) * length(datasets.connected.pre))
df.nm2$sigmas = rep(rep(sigmas, each = length(tried)), times = length(datasets.connected.pre))
df.nm2$dataset = rep(datasets.connected.pre, each = length(sigmas) * length(tried))

tmp = merge(df.nm2, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm2.pvalue = ddply(tmp, .(dataset, sigmas), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / 500, 
                      zvalue.lev = 1- pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)),
                      pvalue.W2 = sum(w2.x > w2.y) / 500, 
                      zvalue.W2 = 1- pnorm((mean(w2.y) - mean(w2.x)) / sd(w2.x)),
                      pvalue.Hq = sum(Hq.x > Hq.y) / 500, 
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
    B = swaplinks(A, HowManyToTry = 2000)  # swaping links that keep the degree distribution constant
    get.measures(B)
  })
})
df.nm3$tried = rep(tried, times = length(datasets.connected.pre))
df.nm3$dataset = rep(datasets.connected.pre, each = length(tried))

df.nm3.random = df.nm3  # backup 
df.nm3 = merge(df.nm3, df.empirical, by = c('dataset'))  # merge two frames by column 'dataset'
df.nm3.pvalue = ddply(df.nm3, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / 500, 
                      zvalue.lev = 1 - pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)))
sum(df.nm3.pvalue$pvalue.lev > 0.95 | df.nm3.pvalue$pvalue.lev < 0.05)  # pvalues not in [0.05, 0.95]

# get the correlation between Hq and lev
colnames(df.nm3)[15] = 'WNODF.x'; colnames(df.nm3)[31] = 'WNODF.y'  # rename the column name of 'weighted NODF' to 'WNODF'
df.nm3.corr = ddply(df.nm3, .(dataset), summarise, pearson.Hq = cor(sqrt(sqrt(Hq.x)), lev.weight.x),  # peason coefficient between Hq and lev
                    R2.Hq = summary(lm(sqrt(sqrt(Hq.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between Hq and lev
                    pearson.NODF = cor(sqrt(sqrt(NODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.NODF = summary(lm(sqrt(sqrt(NODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.WNODF = cor(sqrt(sqrt(WNODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.WNODF = summary(lm(sqrt(sqrt(WNODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.WINE = cor(sqrt(sqrt(wine.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.WINE = summary(lm(sqrt(sqrt(wine.x)) ~ lev.weight.x, data=faithful))$r.squared)  # R^2 coefficient between NODF and lev
                    
df.nm3.corr.long = melt(df.nm3.corr, id.vars = c('dataset'), measure.vars = c('R2.Hq', 'R2.NODF', 'R2.WNODF', 'R2.WINE'),
                      variable.name = 'corr.type', value.name = 'corr.value')


p <- ggplot(data=df.nm3.corr.long, aes(x = corr.type, y = corr.value))
p + geom_boxplot() +
#  facet_wrap(~nullmodels, scale='free') +
  scale_x_discrete('Different Measures') +
  scale_y_continuous(expression(paste('Correlations with  ', lambda[1]))) +
  theme_bw()


########### measures of Null Model 4 #############
df.nm4 = ldply(datasets.connected.pre, .parallel = TRUE, function(A) {
  A = get(A)
  edges = length(A[A > 0])  # get number of links
  nodes = dim(A)  # number of nodes of two groups
  ldply(tried, .parallel = TRUE, function(i) {
    print(paste(i))
    # B <- matrix(sample(A), nodes[1], nodes[2])
    # generate  a ER random bipartite graph with same number of nodes and edges with empirical networks
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
df.nm4.pvalue = ddply(df.nm4, .(dataset), summarise, pvalue.lev = sum(lev.weight.x > lev.weight.y) / 500, 
                      zvalue.lev = 1 - pnorm((mean(lev.weight.y) - mean(lev.weight.x)) / sd(lev.weight.x)))
sum(df.nm4.pvalue$pvalue.lev > 0.95 | df.nm4.pvalue$pvalue.lev < 0.05)  # pvalues not in [0.05, 0.95]

# get the correlation between Hq and lev
colnames(df.nm4)[15] = 'WNODF.x'; colnames(df.nm4)[31] = 'WNODF.y'  # rename the column name of 'weighted NODF' to 'WNODF'
df.nm4.corr = ddply(df.nm4, .(dataset), summarise, pearson.Hq = cor(sqrt(sqrt(Hq.x)), lev.weight.x),  # peason coefficient between Hq and lev
                    R2.Hq = summary(lm(sqrt(sqrt(Hq.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between Hq and lev
                    pearson.NODF = cor(sqrt(sqrt(NODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.NODF = summary(lm(sqrt(sqrt(NODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.WNODF = cor(sqrt(sqrt(WNODF.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.WNODF = summary(lm(sqrt(sqrt(WNODF.x)) ~ lev.weight.x, data=faithful))$r.squared,  # R^2 coefficient between NODF and lev
                    pearson.WINE = cor(sqrt(sqrt(wine.x)), lev.weight.x),  # peason coefficient between NODF and lev
                    R2.WINE = summary(lm(sqrt(sqrt(wine.x)) ~ lev.weight.x, data=faithful))$r.squared)  # R^2 coefficient between NODF and lev

df.nm4.corr.long = melt(df.nm4.corr, id.vars = c('dataset'), measure.vars = c('R2.Hq', 'R2.NODF', 'R2.WNODF', 'R2.WINE'),
                        variable.name = 'corr.type', value.name = 'corr.value')

p <- ggplot(data=df.nm4.corr.long, aes(x = corr.type, y = corr.value))
p + geom_boxplot() +
  #  facet_wrap(~nullmodels, scale='free') +
  scale_x_discrete('Different Measures') +
  scale_y_continuous(expression(paste('Correlations with  ', lambda[1]))) +
  theme_bw()



