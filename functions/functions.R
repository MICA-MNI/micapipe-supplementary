# R script
#
# Functionf for micapipe figures 
# 
# R version 3.6.3
#
# Created by RRC on December 2021 (the second year of the pademic)

# --------------------------------------------------------------------------------- #

# require
require("igraph")

#### Funtions ####
# function: Strenght
strenght_wu <- function(M){
  S <- apply(M,2,sum)
  return(S)
}

# function: Weighted charachteristical path length
nodal_lenght <- function(M){
  M <- 1/M
  M[which(M==-Inf)] <- Inf
  Lnet <- graph_from_adjacency_matrix( (M) ,mode="undirected",weighted = TRUE,diag = FALSE) 
  D <- distances(Lnet,algorithm = "dijkstra",weights=E(Lnet)$weight)
  if(sum(is.infinite(D))>0) D[is.infinite(D)] <- max(D[!is.infinite(D)]) # in order to avoid infinite numbers (Fornito et al., 2010)
  diag(D) <- 0
  L <- apply(D,1,sum)/(dim(M)[1]-1)
  return(L)
}

# function: Weighted clustering coefficent
cc_wu <- function(M){
  Mnet <- graph_from_adjacency_matrix( M ,mode="undirected",weighted = TRUE,diag = FALSE) 
  CCw <- transitivity(Mnet,type="weighted",weights=E(Mnet)$weight)
  return(CCw)
}

# function: fischer transformation
fisher <- function(rho) {
  r = (1+rho)/(1-rho)
  z = 0.5*log(r,base = exp(1))
  return(z)
}

# function: load connectomes
load.mtx <- function(File, conn) {
  # Load the cortical connectome
  mtx <- as.matrix(read.csv(File, sep=" ", header <- FALSE,))
  n <- dim(mtx)[1]
  
  if (conn!='GD') {
    # print("Filling the lower triangle of the matrix")
    mtx[lower.tri(mtx)] <- t(mtx)[lower.tri(mtx)]
  }
  
  if (conn=='SC' | conn=='FC') {
    mtx <- mtx[50:n, 50:n]
    indx <- ((dim(mtx)[1]-1)/2)+1
    mtx <- mtx[-indx,-indx]
  }
  
  if (conn=='MPC' | conn=='GD' ) {
    mtx <- mtx[-1,-1]
    indx <- ((dim(mtx)[1]-1)/2)+1
    mtx <- mtx[-indx,-indx]
    
    if (conn=='FC') {
      print("Fisher transform")
      mtx <- fisher(mtx)
      # replace inf with 0
      mtx[is.infinite(mtx)] <- 0
    }
  }
  return(mtx)
}

# function: Proportional threshold
proportional.thr <- function(M, thr=0.8) {
  M[which(M < quantile(M, prob=thr), arr.ind=TRUE)] <-0
  return(M)
}

# function: Normalize from 0 to 1
norm01 <- function(x) {
  zi <- (x -min(x)) / (max(x)-min(x))
  return(zi)
}
