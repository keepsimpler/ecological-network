library(vegan)
library(bipartite)
library(igraph)

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

estimate.dominant.eigenvalue <- function(A) {
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

order.by.rowsums.and.colsums <- function(comm, decreasing = FALSE) {
  rfill <- rowSums(comm)
  cfill <- colSums(comm)
  rorder <- order(rfill, decreasing = decreasing)
  corder <- order(cfill, decreasing = decreasing)
  comm <- comm[rorder, corder]
  comm
}

###############################################################################
### nestedness introduced by Bastolla and collaborators based on ComMon NeighBors
### input: incident matrix
nest.cmnb <- function (comm)
{
  comm <- ifelse(comm > 0, 1, 0)  # get binary network
  rfill <- rowSums(comm)
  cfill <- colSums(comm)
  if (any(rfill == 0) || any(cfill == 0)) 
    stop('can not have zero rows and cols!')
  nr <- NROW(comm)
  nc <- NCOL(comm)
  fill <- sum(rfill)/prod(dim(comm))  # connectence

  overlapP <- comm %*% t(comm)  # overlap matrix of Plants
  tmp1 = matrix(rep(rfill, times = nr), byrow = T, ncol = nr) 
  tmp2 = matrix(rep(rfill, times = nr), byrow = F, ncol = nr)
  tmp = pmin(tmp1, tmp2)
  overlapP = overlapP / tmp
  overlapA <- t(comm) %*% comm  # overlap matrix of Animals
  tmp1 = matrix(rep(cfill, times = nc), byrow = T, ncol = nc) 
  tmp2 = matrix(rep(cfill, times = nc), byrow = F, ncol = nc)
  tmp = pmin(tmp1, tmp2)
  overlapA = overlapA / tmp

  N.rows <- mean(overlapP[upper.tri(overlapP)])
  N.columns <- mean(overlapA[upper.tri(overlapA)])
  CMNB <- (sum(c(overlapP[upper.tri(overlapP)], overlapA[upper.tri(overlapA)])))/
    ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
  out <- list(fill = fill, N.columns = N.columns, N.rows = N.rows, CMNB = CMNB)
  out
}

nest.cmnb2 <- function (comm)
{
  comm <- ifelse(comm > 0, 1, 0)  # get binary network
  rfill <- rowSums(comm)
  cfill <- colSums(comm)
  if (any(rfill == 0) || any(cfill == 0)) 
    stop('can not have zero rows and cols!')
  nr <- NROW(comm)
  nc <- NCOL(comm)
  fill <- sum(rfill)/prod(dim(comm))  # connectence

  overlapP <- comm %*% t(comm)  # overlap matrix of Plants
  tmp1 = matrix(rep(rfill, times = nr), byrow = T, ncol = nr) 
  tmp2 = matrix(rep(rfill, times = nr), byrow = F, ncol = nr)
  tmp = pmin(tmp1, tmp2)
  N.rows = sum(overlapP[upper.tri(overlapP)]) / sum(tmp[upper.tri(tmp)])
  overlapA <- t(comm) %*% comm  # overlap matrix of Animals
  tmp1 = matrix(rep(cfill, times = nc), byrow = T, ncol = nc) 
  tmp2 = matrix(rep(cfill, times = nc), byrow = F, ncol = nc)
  tmp = pmin(tmp1, tmp2)
  N.columns = sum(overlapA[upper.tri(overlapA)]) / sum(tmp[upper.tri(tmp)])
  
  CMNB <- (N.rows + N.columns)/2
  out <- list(fill = fill, N.columns = N.columns, N.rows = N.rows, CMNB = CMNB)
  out
}

###############################################################################
#### nestedness defination of NODF
#### Almeida-Neto, M., Guimarães, P., Guimarães, P. R., Loyola, R. D. & Ulrich, W.
#### A consistent metric for nestedness analysis in ecological systems: 
#### reconciling concept and measurement. Oikos 117, 1227–1239 (2008).

nest.nodf <- function (comm, order = TRUE)
{
  comm <- ifelse(comm > 0, 1, 0)  # get binary network
  rfill <- rowSums(comm)
  cfill <- colSums(comm)
  if (order) {
    rorder <- order(rfill, decreasing = TRUE)
    corder <- order(cfill, decreasing = TRUE)
    comm <- comm[rorder, corder]
    rfill <- rfill[rorder]
    cfill <- cfill[corder]
  }
  nr <- NROW(comm)
  nc <- NCOL(comm)
  fill <- sum(rfill)/prod(dim(comm))  # connectence
  N.paired.rows <- numeric(nr * (nr - 1)/2)
  N.paired.cols <- numeric(nc * (nc - 1)/2)
  counter <- 0
  for (i in 1:(nr - 1)) {
    first <- comm[i, ]
    for (j in (i + 1):nr) {
      counter <- counter + 1
      if (rfill[i] <= rfill[j] || any(rfill[c(i, j)] == 0))
        next
      N.paired.rows[counter] <- sum(first + comm[j, ] == 2)/rfill[j]
    }
  }
  counter <- 0
  for (i in 1:(nc - 1)) {
    first <- comm[, i]
    for (j in (i + 1):nc) {
      counter <- counter + 1
      if (cfill[i] <= cfill[j] || any(cfill[c(i, j)] == 0))
        next
      N.paired.cols[counter] <- sum(first + comm[, j] == 2)/cfill[j]
    }
  }
  N.columns <- mean(N.paired.cols)
  N.rows <- mean(N.paired.rows)
  NODF <- (sum(c(N.paired.rows, N.paired.cols)))/
    ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
  out <- list(fill = fill, N.columns = N.columns, N.rows = N.rows, NODF = NODF)
  # class(out) <- "nestednodf"
  out
}

## NODF which doesn't depend on the sorting of rows and cols
nest.nodf2 <- function (comm, order = TRUE)
{
  comm <- ifelse(comm > 0, 1, 0)  # get binary network
  rfill <- rowSums(comm)
  cfill <- colSums(comm)
  if (any(rfill == 0) || any(cfill == 0)) 
    stop('can not have zero rows and cols!')
  if (order) {
    rorder <- order(rfill, decreasing = TRUE)
    corder <- order(cfill, decreasing = TRUE)
    comm <- comm[rorder, corder]
    rfill <- rfill[rorder]
    cfill <- cfill[corder]
  }
  nr <- NROW(comm)
  nc <- NCOL(comm)
  fill <- sum(rfill)/prod(dim(comm))  # connectence
  N.paired.rows <- numeric(nr * (nr - 1)/2)
  N.paired.cols <- numeric(nc * (nc - 1)/2)
  counter <- 0
  for (i in 1:(nr - 1)) {
    for (j in (i + 1):nr) {
      counter <- counter + 1
      if (rfill[i] == rfill[j]) next  # why??
      N.paired.rows[counter] <- (comm[i, ] %*% comm[j, ])/min(rfill[i], rfill[j])
    }
  }
  counter <- 0
  for (i in 1:(nc - 1)) {
    for (j in (i + 1):nc) {
      counter <- counter + 1
      if (cfill[i] == cfill[j]) next  # why??
      N.paired.cols[counter] <-  (comm[, i] %*% comm[, j])/min(cfill[i], cfill[j])
    }
  }
  N.columns <- mean(N.paired.cols)
  N.rows <- mean(N.paired.rows)
  NODF <- (sum(c(N.paired.rows, N.paired.cols)))/
    ((nc * (nc - 1)/2) + (nr * (nr - 1)/2))
  out <- list(fill = fill, N.columns = N.columns, N.rows = N.rows, NODF = NODF)
  # class(out) <- "nestednodf"
  out
}
