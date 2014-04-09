##### evaluate 'Disentangling nestedness from models of ecological complex'

out = lv.1.check(dataset = Safariland.C.Pre)
matplot(t, lvout[, 2:10], type = "l", lwd = 1.5, col = 1)

############ Nestedness definition
source('nestedness.R')  # import the functions defining different nestednesss
for (i in datasets.connected.pre) {
  A = get(i)
  nest1 = nodf.old(A, order = F)['NODF']
  nest2 = nestednodf(A,order = F)$statistic['NODF']
  print(paste(nest1, nest2))
}

A = Safariland.C.Pre  # incident matrix
A[A > 0] = 1  # transform to binary matrix
B = as.one.mode(A)  # transform to adjacency matrix
B2 = B %*% B # get the overlap matrix, which indicate the common neighbors between species
sum(B2)  # the sum of all the common neighbors between species
degrees = rowSums(B)  # the degrees (number of neighbors) of species
degrees^2  # the number of neighbors of neighbors of species
print(paste(sum(B2), sum(degrees^2)))  # 
