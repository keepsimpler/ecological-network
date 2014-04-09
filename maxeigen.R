est.max.eigen <- function(A) {
  A = GetAdjacencyMatrix(A)
  B = A
  B[B!=0] = 1
  eigenvalues = eigen(B)$values
  max.eigenvalues = eigenvalues[1]
  BM2 =  sum(diag(B %*% B))
  #BM3 =  sum(diag(B %*% B %*% B)) / s
  BM4 =  sum(diag(B %*% B %*% B %*% B))
  s = dim(B)[1]
  lambda1 = sqrt( ( (s-2)*BM4 + 8 * BM2^2 / (s+2) - 2*BM2^2 ) / (2*(s+2))) + 2*BM2 / (s+2)
  lambda1 = sqrt(lambda1)
  
  r = (BM2 - 2*lambda1^2) / (s - 2)
  
  strengths2 = rowSums(B^2)
  q2 = sum(strengths2^2)
  
  chung = sum(rowSums(B)**2) / sum(rowSums(B))
  
  c(max.eigenvalues = max.eigenvalues, lambda1 = lambda1, chung = chung)
}

## get the adjacency matrix according to the interaction frequency matrix
GetAdjacencyMatrix <- function(A){
  NumP <- dim(A)[1]
  NumA <- dim(A)[2]
  S <- NumP + NumA
  Adj <- matrix(0, S, S)
  Adj[1:NumP, (NumP + 1):S] <- A
  Adj <- Adj + t(Adj)
  return(Adj)
}
