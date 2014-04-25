library(vegan)
library(bipartite)
library(igraph)

###############################################################################
#' @title Generate a connected graph using package [igraph]
#'
#' @param s size of network. 
#' if graph type is bipartite, s[1], s[2] represent size of two groups; else s is size of network
#' @param k average degree for the network.
#' 1 < k < s for unipartite network, 1 < k < s[1]*s[2]/(s[1]+s[2]) for bipartite network.
#' @param gtype Graph type generated: 'bipartite'
#' @param ... the parms conform to the implementation functions of [igraph]
#' @return the connected graph
#' @details .  
#' @import igraph
graph.connected <- function(s, k, gtype, ...) {
  library(igraph)
  if (gtype == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO node groups should be designated. we have assumed second size to be the first size.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  repeat {  # generate a connected graph
    if (gtype == 'bipartite') {
      G = bipartite.random.game(s[1], s[2], type = 'gnm', m = k * (s[1] + s[2]))
    }
    else {
      G = erdos.renyi.game(s, m = k * s, type = 'gnm')
    }
    if (igraph::is.connected(G)) break  # until a connected graph is generated
  }
  G
}

###############################################################################
#' @title Nestedness optimization algrithm by rewiring links to nodes with more links. (richer get richer)
#'
#' @param B incidence matrix of bipartite network, rows and cols represent two groups of nodes/species
#' @param HowManyToTry the times to try for rewiring links
#' @return the incidence matrix whose nestedness has been optimized by rewiring links.
#' @details .  
#' @import bipartite
rewirelinks.richer <- function(B, HowManyToTry) {
  library(bipartite)
  B = sortweb(B)  # sort rows and cols descending, ensure the chosen species later has more interactions
  count1 <- 0
  NumP <- dim(B)[1]
  NumA <- dim(B)[2]
  while (count1 < HowManyToTry){
    count1 <- count1 + 1
    ## pick one interaction between two random species
    repeat {
      row1 <- sample(1:NumP, 1)
      col1 <- sample(1:NumA, 1)
      if (B[row1, col1] != 0) break
    }
    ## random choose another species
    if (runif(1) < 0.5) {  # choose another plant
      row2 =  sample(1:row1, 1)  # choose random plant with more interactions which is entured by [sortweb]
      # Three exceptions: 1. the new chosen species [row2] is same with the old species [row1]
      # 2. the new chosen species [row2] already has interaction with [col1]
      # 3. the old species [row1] has only one interaction.
      if (row2 < row1 && B[row2, col1] == 0 && sum(B[row1,]) > 1) {
        B[row2, col1] = B[row1, col1]
        B[row1, col1] = 0
      }
    }
    else {  # choose another animal
      col2 =  sample(1:col1, 1)
      if (col2 < col1 && B[row1, col2] == 0 && sum(B[,col1]) > 1 ) {
        B[row1, col2] = B[row1, col1]
        B[row1, col1] = 0
      }
    }
    sortweb(B)  # sort rows and cols descending for the next link rewiring.
  }
  B
}

###############################################################################
#' @title Swap links of bipartite network, that will keep the node degree distribution.
#'
#' @param B incidence matrix of bipartite network, rows and cols represent two groups of nodes/species
#' @param HowManyToTry the times to try for swapping links
#' @return the incidence matrix whose links being randomly swapped.
#' @details .  
swaplinks <- function(B, HowManyToTry = 5000) {
  count1 <- 0
  NumP <- dim(B)[1]
  NumA <- dim(B)[2]
  while (count1 < HowManyToTry){
    count1 <- count1 + 1
    ## pick two rows and two columns
    row1 <- sample(1:NumP, 1)
    row2 <- sample(1:NumP, 1)
    col1 <- sample(1:NumA, 1)
    col2 <- sample(1:NumA, 1)
    ## check swappable
    if (B[row1, col1] == 0.0 && B[row1, col2] > 0.0 && B[row2, col1] > 0.0 && B[row2, col2] == 0.0){
      ## swap
      B[row1, col1] <- B[row1, col2]
      B[row1, col2] <- 0.0
      B[row2, col2] <- B[row2, col1]
      B[row2, col1] <- 0.0
    }
    else{
      if (B[row1, col1] > 0.0 && B[row1, col2] == 0.0 && B[row2, col1] == 0.0 && B[row2, col2] > 0.0){
        ## swap
        B[row1, col2] <- B[row1, col1]
        B[row1, col1] <- 0.0
        B[row2, col1] <- B[row2, col2]
        B[row2, col2] <- 0.0
      }
    }
  }
  return(B)
}


# order incidence matrix of bipartite network by ascending or deascending of rows and cols totals.
# Replaced by bipartite::sortweb
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
