# R script
#
# Figure 4 micapipe 
# Datasets-similarities, correlation of the mean edges, and graph theoretical analisys
# R version 3.6.3
#
# Created by RRC on December 2021 (the second year of the pademic)

# Set the environment
require("RColorBrewer")
require("viridis")
library("igraph")
library("tidyverse")
library("corrplot")
library("fsbrain")
library("freesurferformats")
library("rgl")

# heper functions
source("~/git_here/micapipe-supplementary/functions/functions.R")

# --------------------------------------------------------------------------------- #
####  Set the working directory to your subjec's directory ####
P <- "/Users/rcruces/tmp/micaConn/"
setwd(P)

# List the datasets
datasets <- c('MICS', 'CamCAN', 'EpiC', 'EpiC2', 'audiopath', 'sudmex', 'MSC')
modalities <- c('GD', 'SC', 'FC', 'MPC')

# Empty list of lists of arrays
dataAVG <- list(GD=list(n100=array(0,dim = c(100,100,7)), n400=array(0,dim = c(400,400,7)), n1000=array(0,dim = c(1000,1000,7))),
                SC=list(n100=array(0,dim = c(100,100,7)), n400=array(0,dim = c(400,400,7)), n1000=array(0,dim = c(1000,1000,7))),
                FC=list(n100=array(0,dim = c(100,100,7)), n400=array(0,dim = c(400,400,7)), n1000=array(0,dim = c(1000,1000,7))),
                MPC=list(n100=array(0,dim = c(100,100,7)), n400=array(0,dim = c(400,400,7)), n1000=array(0,dim = c(1000,1000,7))))

# --------------------------------------------------------------------------------- #
#### Read all the subjects #### 
# For each dataset and modality
for (k in 1:length(datasets)){
  dataset <- datasets[k]
  
  start.time <- Sys.time()
  
  for (i in 1:length(modalities)) {
    mod <- modalities[i]
    
    # list files
    Files <- list.files(path = paste0(P, mod, '/', dataset), pattern = ".txt")
    
    # Only runs if connectomes are found
    if (length(Files)!=0) {
      print(paste('\n RUNNING dataset', dataset,mod))
      
      # for each parcellation
      parcs <- c('100', '400', '1000')
      
      for (j in 1:length(parcs)) {
        parc <- parcs[j]
        
        # load each parcellation into a separate array
        files.parc <- Files[grep(paste0('schaefer-',parc,'_'), Files)]
        
        # Dimensions of array
        n <- dim(load.mtx(paste0(mod,"/",dataset,'/',files.parc[1]), conn=mod ))[1]
        N <- length(files.parc)
        
        # Build array of all matrices
        print(paste0("[INFO] ... reading parcellation: ", parc, ", dataset: ", dataset, "-",mod))
        mtxs <- array(0,dim = c(n,n,N))
        for (w in 1:N) {
          mtxs[,,w] <- load.mtx(paste0(mod,"/",dataset,'/',files.parc[w]), conn = mod)
        }
        
        # Calculate the mean
        Mavg <- apply(mtxs,c(1,2), mean)
        
        # Assign matrix to object
        dataAVG[[i]][[j]][,,k] <- Mavg
        
        # Computational Time
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
      }
      
    } else {
      print(paste('[INFO].... NOTHING found connectomes in', dataset,mod))
    }
    
  }
}

# save(dataAVG, file = 'all_modalities-avg.RData')

# --------------------------------------------------------------------------------- #
#### Datasets correlations #### 
load(file = 'all_modalities-avg.RData')

# Granularity (only 400 for the figure 4)
grain <- 400

# ------------------------------------------------------------------------------- #
# Function: plot correlation matrix from dataframe
plot.corr <- function(Data, Title) {
  corr <- cor(Data); corr[is.na(corr)] <- 0; corr <- abs(corr)
  svg(paste0(Title, "_datasets-correlation.svg"), width = 4.5, height = 4.5)
  g <- corrplot(corr, is.corr = TRUE, method = 'color', type = 'lower', diag = FALSE, title = Title, 
           cl.lim = c(0,1),cl.ratio = 0.6, cl.pos = 'b', col = rep(brewer.pal(9,'Reds'),2), tl.col = 'black', addgrid.col = NA)
  dev.off()
  return(g)
}

# Between datasets correlation matrices
for (mod in modalities) {
  # N databases
  M <- dataAVG[[mod]][[paste0('n',grain)]]
  N <- dim(M)[3]
  n <- dim(M)[1]
  
  # ------------------------------------ #
  # Edges
  E <- sapply(1:N, function(x) c(M[,,x][upper.tri(M[,,x])]) )
  E <- data.frame(E); colnames(E) <- datasets
  print(paste0('[INFO].... Edges: ',mod,'-',grain))
  plot.corr(E,  paste0('edges_',mod,'-',grain))
  
  # ------------------------------------ #
  # Strength
  # Sparcity threshold, top 20% of the connections except in SC
  if (mod!='SC') {
    tmp <- sapply(1:N, function(x) proportional.thr(M[,,x], 0.8))
    dim(tmp) <- dim(M)
    M <- tmp; rm(tmp)
  }
  
  # Remove negative values for FC
  if (mod=='FC') { M[which(M<0)] <- 0 }
  if (mod=='MPC') { M <- abs(M) }
  
  # plot
  S <- sapply(1:N, function(x) strenght_wu(M[,,x]))
  S <- data.frame(S); colnames(S) <- datasets
  print(paste0('[INFO].... Strength: ',mod,'-',grain))
  plot.corr(S, paste0('strength_',mod,'-',grain))
  
  # ------------------------------------ #
  # Characteristic path length
  if (mod=='GD') {
    # left hemisphere
    L.L <- sapply(1:N, function(x) nodal_lenght(M[1:(n/2),1:(n/2),x]))
    # right hemisphere
    L.R <- sapply(1:N, function(x) nodal_lenght(M[((n/2)+1):n,((n/2)+1):n,x]))
    L <- rbind(L.L, L.R)
  } else {
    L <- sapply(1:N, function(x) nodal_lenght(M[,,x]))
  }
  L <- data.frame(L); colnames(L) <- datasets
  print(paste0('[INFO].... Path length: ',mod,'-',grain))
  plot.corr(L, paste0('path-length_',mod,'-',grain))
  
  # ------------------------------------ #
  # Cluster coefficient
  if (mod=='SC') {
    tmp <- sapply(1:N, function(x) proportional.thr(M[,,x], 0.8))
    dim(tmp) <- dim(M)
    M <- tmp; rm(tmp)
  }
  if (mod=='GD') {
    # left hemisphere
    CCw.L <- sapply(1:N, function(x) cc_wu(M[1:(n/2),1:(n/2), x]))
    CCw.L[which(is.na(CCw.L))] <- 0
    # right hemisphere
    CCw.R <- sapply(1:N, function(x) cc_wu(M[((n/2)+1):n,((n/2)+1):n, x]))
    CCw.R[which(is.na(CCw.R))] <- 0
    CCw <- rbind(CCw.L, CCw.R)
  } else {
    CCw <- sapply(1:N, function(x) cc_wu(M[,,x]))
    CCw[which(is.na(CCw))] <- 0
  }
  CCw <- data.frame(CCw); colnames(CCw) <- datasets
  print(paste0('[INFO].... Cluster Coeff: ',mod,'-',grain))
  plot.corr(CCw, paste0('cluster-coeff_',mod,'-',grain))
  
}
