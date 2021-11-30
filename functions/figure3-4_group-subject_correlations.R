# R script
#
# Figure 3 and 4 micapipe 
# Subject-Group mean Correlation coefficients, edges, and graph theoretical analisys
# R version 3.6.3
#
# Created by RRC on September 2021 (the second year of the pademic)

# Set the environment
require("RColorBrewer")
require("viridis")
library("igraph")
library("tidyverse")

# --------------------------------------------------------------------------------- #
#### Definen functions ####
# Strenght
strenght_wu <- function(M){
  S <- apply(M,2,sum)
  return(S)
}

# Weighted charachteristical path length
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

# Weighted clustering coefficent
cc_wu <- function(M){
  Mnet <- graph_from_adjacency_matrix( M ,mode="undirected",weighted = TRUE,diag = FALSE) 
  CCw <- transitivity(Mnet,type="weighted",weights=E(Mnet)$weight)
  return(CCw)
}

#----------------------------------------------------------#
# Coeficiente de Clustering - BINARY UNDIRECTED
clustering_coef_bu <- function(A){
  A[A > 0] = 1 # Make sure is a binary matrix
  n <- dim(A)[1]
  C <- rep(0,n) # QUe hacer en caso de que no tenga CC?
  
  for (u in 1:n) {
    Ver <- which(A[u,]!=0)
    k <- length(Ver)
    # degree must be at least 2
    if (k >= 2) {
      S <- A[Ver,Ver]
      C[u] = sum(S[,])/(k^2-k)
    }                 
  }
  return(C)
}

fisher <- function(rho) {
  r = (1+rho)/(1-rho)
  z = 0.5*log(r,base = exp(1))
  return(z)
}

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

proportional.thr <- function(M, thr=0.8) {
  M[which(M < quantile(M, prob=thr), arr.ind=TRUE)] <-0
  return(M)
}

# --------------------------------------------------------------------------------- #
#### Group - subject correlations ####
# Set the working directory to your subjec's directory
setwd("/Users/rcruces/tmp/micaConn/")

# List the datasets
datasets <- c('MICS', 'CamCAN', 'EpiC', 'EpiC2', 'audiopath', 'sudmex', 'MSC')
modalities <- c('GD', 'SC', 'FC', 'MPC')

# For each dataset and modality
for (dataset in datasets){
  
  # Create empty variables for each dataset
  corr <- c()
  means <- c()
  
  for (mod in modalities) {
    
    # list files
    Files <- list.files(path = paste0("/Users/rcruces/tmp/micaConn/", mod, '/', dataset), pattern = ".txt")
    
    # Only runs if connectomes are found
    if (length(Files)!=0) {
      print(paste('[INFO].... RUNNING dataset', dataset,mod))
      
      # for each parcelation
      # parcs <- c('100', '400', '1000')
      parcs <- c('100', '400')
      
      for (parc in parcs) {
        # load each parcellation into a separate array
        files.parc <- Files[grep(paste0('schaefer-',parc,'_'), Files)]
        
        # Dimensions of array
        n <- dim(load.mtx(paste0(mod,"/",dataset,'/',files.parc[1]), conn=mod ))[1]
        N <- length(files.parc)
        
        # Build array of all matrices
        print(paste0("[INFO] ... reading parcellation: ", parc, ", dataset: ", dataset, "-",mod))
        mtxs <- array(0,dim = c(n,n,N))
        for (i in 1:N) {
          mtxs[,,i] <- load.mtx(paste0(mod,"/",dataset,'/',files.parc[i]), conn = mod)
        }
        
        # Calculate the mean
        Mavg <- apply(mtxs,c(1,2), mean)
        
        # Plot the log matrix
        if (mod=='SC') {
          Col=brewer.pal(9,'Purples')
          Mplot <- log(Mavg)
        } else if (mod=='FC') {
          Col=brewer.pal(9,'Reds')
          Mplot <- Mavg
        } else if (mod=='GD') {
          Col=brewer.pal(9,'Blues')
          Mplot <- Mavg    
        } else if (mod=='MPC') {
          Col=brewer.pal(9,'Greens')
          Mplot <- Mavg    
        }
        png(paste0(mod, "_", dataset, "_", parc, ".png"), units = "px", width = 500, height = 500)
        image(Mplot, axes=FALSE, main=paste(mod, "mean, granularity:", parc), col=Col)
        dev.off()
        
        # ------------------------------------------------------------------------------- #
        #### EDGES: Subject to group correlations #### 
        # Correlation of the edges | upper triangle only
        Cor <- c()
        for (i in 1:N) {
          Cor <- c(cor(c(mtxs[,,i][upper.tri(mtxs[,,i])]), c(Mavg[upper.tri(Mavg)])), Cor)
        }
        
        # ------------------------------------------------------------------------------- #
        #### Graph theoretical analysis: Subject to group correlations ####
        # Sparcity threshold, top 20% of the connections except in SC
        if (mod!='SC') {
          tmp <- sapply(1:N, function(x) proportional.thr(mtxs[,,x], 0.8))
          dim(tmp) <- dim(mtxs)
          mtxs <- tmp; rm(tmp)
        }
        
        # Remove negative values for FC
        if (mod=='FC') { mtxs[which(mtxs<0)] <- 0 }
        
        # Strenght
        print("[INFO] ... calculating Nodal strength")
        S <- t(sapply(1:N, function(x) strenght_wu(mtxs[,,x])))
        
        # Mean strenght
        Savg <- apply(S,2,mean)
        
        # Correlation of the strenght
        Scor <- c()
        for (i in 1:N) {
          Scor <- c(cor(Savg, S[i,]), Scor)
        }
        
        # ------------------------------------------------------------------------------- #
        # Cluster coefficient
        print("[INFO] ... calculating Weighted Clustering Coefficient")
        if (mod=='GD') {
          # left hemisphere
          CCw.L <- t(sapply(1:N, function(x) cc_wu(mtxs[1:(n/2),1:(n/2), x])))
          CCw.L[which(is.na(CCw.L))] <- 0
          # right hemisphere
          CCw.R <- t(sapply(1:N, function(x) cc_wu(mtxs[((n/2)+1):n,((n/2)+1):n, x])))
          CCw.R[which(is.na(CCw.R))] <- 0
          CCw <- cbind(CCw.L, CCw.R)
        } else {
          CCw <- t(sapply(1:N, function(x) cc_wu(mtxs[,,x])))
          CCw[which(is.na(CCw))] <- 0
        }
        
        # Mean CCw
        Cavg <- apply(CCw,2,mean)
        
        # Correlation of the cluster coefficient
        Ccor <- c()
        for (i in 1:N) {
          Ccor <- c(cor(Cavg, CCw[i,]), Ccor)
        }
        
        # ------------------------------------------------------------------------------- #
        # Path
        if (mod=='MPC') { mtxs <- abs(mtxs)}
        print("[INFO] ... calculating Path lenght")
        if (mod=='GD') {
          # left hemisphere
          L.L <- t(sapply(1:N, function(x) nodal_lenght(mtxs[1:(n/2),1:(n/2),x])))
          # right hemisphere
          L.R <- t(sapply(1:N, function(x) nodal_lenght(mtxs[((n/2)+1):n,((n/2)+1):n,x])))
          L <- cbind(L.L, L.R)
        } else {
          L <- t(sapply(1:N, function(x) nodal_lenght(mtxs[,,x])))
        }
        
        # Mean P
        Lavg <- apply(L,2,mean)
        
        # Correlation of the path length
        Lcor <- c()
        for (i in 1:N) {
          Lcor <- c(cor(Lavg, L[i,]), Lcor)
        }
        
        # Create a data.frame with the results
        corr <- rbind(corr, data.frame(granularity=rep(parc,N), modality=rep(mod,N), edges=Cor, strength=Scor, path=Lcor, cluster=Ccor))
        means <- rbind(means, data.frame(roi=1:n, granularity=rep(parc,n), modality=rep(mod,n), strength=Savg, path=Lavg, cluster=Cavg))
      }
      
    } else {
      print(paste('[INFO].... NOTHING found connectomes in', dataset,mod))
    }
    
  }
  
  # save the Modality by dataset
  write.csv(corr, file=paste0(dataset, "_correlations.csv"), row.names = FALSE)
  write.csv(means, file=paste0(dataset, "_node-means.csv"), row.names = FALSE)
}

