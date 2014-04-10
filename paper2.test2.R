


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


