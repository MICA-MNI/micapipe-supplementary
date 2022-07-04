# ----------------------------------------------------------------------
# Iter-Itra subject variability test
# 
# ----------------------------------------------------------------------
# Libraries
source("/host/yeatman/local_raid/rcruces/git_here/micapipe-supplementary/functions/functions.R")
library(corrplot)
library(RColorBrewer)
library(scales)

# ------------------------------------------------ #
#### Intra-subject vs inter-subject variability ####
# Calculate subject variability between two given matrices
# see: Seguin, Caio, Robert E. Smith, and Andrew Zalesky. "Connectome Spatial Smoothing (CSS): concepts, methods, and evaluation." NeuroImage (2022): 118930.
sv <- function(M1, M2){
  # concatenate matrices for faster computation
  M <- cbind(M1, M2)
  Ns <- dim(M)
  # Subject variability
  SV <- mean(apply(M, 1, function(x) cor(x[1:Ns[1]], x[(Ns[1]+1):Ns[2]])), na.rm = TRUE)
  return(SV)
}
# Functions: Identifiability absolute effect size with pooled SD
iden <- function(intra, inter){
  m1 <- intra
  m2 <- inter
  # N size -1
  n1 <- length(m1)-1
  n2 <- length(m2)-1
  # mean difference
  md <- mean(m1)-mean(m2)
  # pooled variance
  S <- sqrt(( n1*var(m1)+n2*var(m2) )/( n1+n2 ))
  D <- abs(md)/S
  # replace Na with 0
  D[is.na(D)] <- 0
  
  return(D)  
}

# ----------------------------------------------------------------------
# GD
# "sub-*/ses-02/anat/surfaces/geo_dist/sub-*_space-fsnative_atlas-schaefer-", schaf,"_GD.txt"
# Structural connnectomes
# "sub-*/ses-02/dwi/connectomes/sub-*_space-dwi_atlas-schaefer-", schaf,"_desc-iFOD2-40M-SIFT2_cor-connectome.txt"
# MPC
# "sub-*/ses-02/anat/surfaces/micro_profiles/sub-*_space-fsnative_atlas-schaefer-", schaf,"_desc-MPC.txt"
# Functional connectomes
# "sub-*/ses-01/func/surfaces/sub-*_rsfmri_space-fsnative_atlas-schaefer-", schaf,"_desc-FC.txt"

# Dataframe with all the paths
df <- data.frame(mod=c("GD", "SC", "MPC", "FC"), 
p1=c("anat/surfaces/geo_dist/sub-*_space-fsnative_atlas-schaefer-","dwi/connectomes/sub-*_space-dwi_atlas-schaefer-","anat/surfaces/micro_profiles/sub-*_space-fsnative_atlas-schaefer-","func/surfaces/sub-*_rsfmri_space-fsnative_atlas-schaefer-"),
p2=c("_GD.txt", "_desc-iFOD2-40M-SIFT2_cor-connectome.txt", "_desc-MPC.txt", "_desc-FC.txt"))
schaf <- '100'

for (k in 4){
  mod=df$mod[k]
  
  # Load files
  #setwd("/data_/mica3/BIDS_EpiC2/derivatives/micapipe")
  setwd("/data_/mica3/BIDS_MICs/derivatives/micapipe")
  files2 <- Sys.glob(paste0("sub-*/ses-02/",df$p1[k], schaf,df$p2[k]))
  subs2 <- gsub("sub-", "", sapply(strsplit(files2, split='/', fixed=TRUE), function(x) (x[1])))
  
  # Dimensions of array
  n <- dim(load.mtx(files2[1], conn=mod ))[1]
  N <- length(files2)
  
  # Build array of all matrices
  mtxs <- array(0,dim = c(n,n,N))
  for (i in 1:N) {
    mtxs[,,i] <- load.mtx(files2[i], conn = mod)
  }
  
  # Load data2
  #setwd("/data_/mica3/BIDS_EpiC/derivatives/micapipe")
  files <- Sys.glob(paste0("sub-*/ses-01/",df$p1[k], schaf,df$p2[k]))
  subs <- gsub("sub-", "", sapply(strsplit(files, split='/', fixed=TRUE), function(x) (x[1])))
  
  # match subjects
  indx <- match(subs2, subs); indx <- indx[!is.na(indx)]
  files <- files[indx]
  subs <- gsub("sub-", "", sapply(strsplit(files, split='/', fixed=TRUE), function(x) (x[1])))
  
  # match missing data from dataset1
  indx2 <- subs2 %in% subs
  subs2 <- subs2[indx2]
  mtxs <- mtxs[,,indx2]
  mean(subs2 %in% subs)
  
  # Build array of all matrices
  N <- dim(mtxs)[3]
  mtxs2 <- array(0,dim = c(n,n,N))
  for (i in 1:N) {
    mtxs2[,,i] <- load.mtx(files[i], conn = mod)
  }
  
  # --------------------------------------------------------------------
  # Create an empty array of subjects x subjecs. Itra-Inter variability
  IIV <- array(0,dim = c(N, N))
  
  # Iterate over each element
  for (i in 1:N){
    for (j in i:N){
      IIV[j,i] <- IIV[i,j] <- sv(mtxs[,,i], mtxs2[,,j])
    }
  }
  png(file=paste0("~/Desktop/iiv_schaefer-",schaf,"_",mod,".png"), width=600, height=350)
  par(mfrow=c(1,2))
  colmin <- 0
  IIV[is.na(IIV)] <- 0
  if (min(IIV)<0){colmin <- min(IIV)}
  corrplot(IIV, method = "color", col = brewer.pal(n=9, name="Reds"), is.corr = FALSE, col.lim = c(colmin, 1), tl.cex = 1, tl.col='black')
  
  hist(IIV[lower.tri(IIV)], breaks =30 , col = "cadetblue2", border = "cadetblue2", main = paste(df$mod[k], schaf), xlab="similarity", probability = FALSE)
  hist(diag(IIV), add=TRUE, col="chocolate1", border="chocolate1", probability = FALSE)
  abline(v=mean(IIV[lower.tri(IIV)]), col="cadetblue4")
  abline(v=mean(diag(IIV)), col="chocolate")
  
  # Reliability (intra-subject)
  Rel <- round(mean(diag(IIV)), 3)
  
  # Uniformity (inter-subjet)
  Uni <- round(mean(IIV[lower.tri(IIV)]), 3)
  
  # Identifiability
  Ide <- round(iden(diag(IIV), IIV[lower.tri(IIV)]), 3)
  
  legend("right", legend=c(paste0("R=",Rel), paste0("U=",Uni), paste0("I=",Ide)), pch=19,col=NA, cex=0.8, bty = "n", )
  dev.off()
}
