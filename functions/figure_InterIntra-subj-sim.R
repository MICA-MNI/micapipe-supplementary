# ----------------------------------------------------------------------
# Iter-Itra subject variability test
# 
# ----------------------------------------------------------------------
# Libraries
source("/host/yeatman/local_raid/rcruces/git_here/micapipe-supplementary/functions/functions.R")
library(corrplot)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)

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
# "sub-*/ses-01/func/surfaces/sub-*_func_space-fsnative_atlas-schaefer-", schaf,"_desc-FC.txt"

# Dataframe with all the paths
df <- data.frame(mod=c("GD", "SC", "MPC", "FC"), 
p1=c("anat/surfaces/geo_dist/sub-*_space-fsnative_atlas-schaefer-","dwi/connectomes/sub-*_space-dwi_atlas-schaefer-","anat/surfaces/micro_profiles/sub-*_space-fsnative_atlas-schaefer-","func/desc-se_task-rest_dir-LR_run-*_bold/surfaces/sub-*_func_space-fsnative_atlas-schaefer-"),
p2=c("_GD.txt", "_desc-iFOD2-20M-SIFT2_cor-connectome.txt", "_desc-MPC.txt", "_desc-FC.txt"))

# Table with values
val <- c()

# Similarity matrices
Smtx <- list()

# Run the Similarity test for each atlas and modality
for (schaf in c('100','400','1000')){
  for (k in 1:4){
    mod=df$mod[k]
    print(paste("Running", mod, schaf))
    # Load files
    setwd("/data_/mica3/BIDS_HCP/test/micapipe")
    files2 <- Sys.glob(paste0("sub-*/",df$p1[k], schaf,df$p2[k]))
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
    setwd("/data_/mica3/BIDS_HCP/retest/micapipe")
    files <- Sys.glob(paste0("sub-*/",df$p1[k], schaf,df$p2[k]))
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

    # Save matrix on list
    Smtx[[paste0(mod,"-",schaf)]] <- IIV
    
    # Create figures
    png(file=paste0("~/Desktop/iiv_schaefer-",schaf,"_",mod,".png"), width=1200, height=700)
    par(mfrow=c(1,2))
    colmin <- 0
    IIV[is.na(IIV)] <- 0
    if (min(IIV)<0){colmin <- min(IIV)}
    if (mod == "GD" ){
      corrplot(IIV, method = "color", col = magma(256), is.corr = FALSE, tl.cex = 1, tl.col='black')
      } else {
    corrplot(IIV, method = "color", col = magma(256), col.lim = c(colmin,1), is.corr = FALSE, tl.cex = 1, tl.col='black')
     }
    
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
    
    # Add legend
    legend("right", legend=c(paste0("R=",Rel), paste0("U=",Uni), paste0("I=",Ide)), pch=19,col=NA, cex=0.8, bty = "n", )
    dev.off()
    
    # Save data to table
    val <- rbind(val, c(schaf, as.character(mod), Rel, Uni, Ide))
  }
}

# Rename the columns of the dataframe with values
colnames(val) <- c("Granularity", "Modality", "Rel", "Uni", "Ide")

# BARPLOT: Intra-subject reliability
barplot(diag(Smtx$`FC-100`), names.arg = subs, las=2, main=paste0("Test-Retest Similarity ",mod,"-",schaf), border=NA, col=ifelse(diag(IIV)<0.6,"gray35","gray65"), ylim = c(0,0.9))
abline(h=0.6, col="red", lty=2)

# MULTI DENSITY PLOT
# Density colors
Col <- c("firebrick3", "darkorange1", "goldenrod")

for (i in 1:4){
  Data <- data.frame(Granularity = factor(rep(c("s100","s400","s1000"), each=53), levels=c("s100","s400","s1000")),
                     Sym = c(diag(Smtx[[i]]), diag(Smtx[[i+4]]), diag(Smtx[[i+8]])))
  Uata <- data.frame(Granularity = factor(rep(c("s100","s400","s1000"), each=1378), levels=c("s100","s400","s1000")),
                     Sym = c(Smtx[[i]][lower.tri(Smtx[[i]])], Smtx[[i+4]][lower.tri(Smtx[[i+4]])], Smtx[[i+8]][lower.tri(Smtx[[i+8]])] ))

  gg <- grid.arrange(
    # Reliability (intra-subject similarity)
    ggplot(data=Data, aes(x=Sym, group=Granularity, fill=Granularity, color=Granularity)) +
    geom_density(adjust=1.5, alpha=.3) +
    scale_fill_manual(values=Col) +
    scale_color_manual(values = Col) +
    scale_x_continuous(limits = c(0, 1)) +
    xlab("Reliability") +
    ylab("Density") +
    ggtitle(df$mod[i]) +
    theme(plot.title = element_text(hjust = 0.5)),
    
    # Uniformity (inter-subject similarity)
    ggplot(data=Uata, aes(x=Sym, group=Granularity, fill=Granularity, color=Granularity)) +
      geom_density(adjust=1.5, alpha=.3) +
      scale_fill_manual(values=Col) +
      scale_color_manual(values = Col) +
      scale_x_continuous(limits = c(0, 1)) +
      xlab("Uniformity") +
      ylab("Density") +
      ggtitle(df$mod[i]) +
      theme(plot.title = element_text(hjust = 0.5)),
    ncol=1, nrow = 2)
  print(gg)
}

# Print all connectome matrices
for ( i in 1:length(names(Smtx))) {
  png(file=paste0("~/Desktop/iiv_",names(Smtx)[i],".png"), width=700, height=700)
  colmin <- 0
  Smtx[[i]][is.na(Smtx[[i]])] <- 0
  if (min(Smtx[[i]])<0){colmin <- min(Smtx[[i]])}
  if (grepl("GD", names(Smtx)[i])){ colmin <- 0.98}
  if (grepl("SC", names(Smtx)[i])){ colmin <- 0.3}
  corrplot(Smtx[[i]], method = "color", col = magma(40), is.corr = FALSE, tl.cex = 1, tl.col='white', main=names(Smtx)[i], col.lim = c(colmin, 1))
  dev.off()
}
val <- data.frame(val)
val[,3:5] <- apply(val[,3:5],2, as.numeric)
val <- data.frame(val)
val$Modality <- factor(val$Modality, levels=c("GD", "SC", "FC", "MPC"))
val$Granularity <- factor(val$Granularity, levels=c("100","400","1000"))

connected_scatter <- function (Y, Title, Ylab="Mean similarity", Print=TRUE, Ylim=c(0,1)) {
  ggp <- ggplot(data = val, aes(x=Granularity, y=Y, group=Modality)) +
    geom_line( aes(color=Modality)) +
    geom_point(shape=21, color="black", size=4, aes(fill=Modality)) +
    ggtitle(Title) +
    ylab(Ylab) +
    scale_y_continuous(limits = Ylim) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15),
          axis.text = element_text(size = 20)) +
    scale_fill_manual(values=c("dodgerblue", "purple","firebrick2", "forestgreen")) +
    scale_color_manual(values = c("dodgerblue", "purple","firebrick2", "forestgreen"))
  if (Print==TRUE){
    print(ggp)
  }
}
grid.arrange(
  connected_scatter(val$Rel, "Reliability", Ylim=c(0.25, 1)),
  connected_scatter(val$Uni, "Uniformity", Ylim=c(0.25, 1)),
  connected_scatter(val$Ide, "Identifiability", Ylab = "Effect size intra-inter", Ylim=c(0,20)),
ncol=3, nrow = 1)