# R script
#
# Figure 3 micapipe 
# Subject-Group mean Correlation coefficients, edges, and graph theoretical analisys
# R version 3.6.3
#
# Created by RRC on December 2021 (the second year of the pademic)

# Set the environment
require("RColorBrewer")
require("viridis")
library("igraph")
library("tidyverse")
library("gridExtra")   # Array for ggplots, equivalent of par(mfrow) 
library("corrplot")
library("fsbrain")
library("freesurferformats")
library("rgl")

# heper functions
source("~/git_here/micapipe-supplementary/functions/functions.R")

# --------------------------------------------------------------------------------- #
#### Group - subject correlations ####
# Set the working directory to your subjec's directory
P <- "/Users/rcruces/tmp/micaConn/"
setwd(P)

# List the datasets
datasets <- c('MICS', 'CamCAN', 'EpiC', 'EpiC2', 'audiopath', 'sudmex', 'MSC')
modalities <- c('GD', 'SC', 'FC', 'MPC')

# For each dataset and modality
for (dataset in datasets){
  # Create empty variables for each dataset
  corr <- c()
  means <- c()
  
  for (mod in modalities) {
    start.time <- Sys.time()
    # list files
    Files <- list.files(path = paste0(P, mod, '/', dataset), pattern = ".txt")
    
    # Only runs if connectomes are found
    if (length(Files)!=0) {
      print(paste('\n RUNNING dataset', dataset,mod))
      
      # for each parcellation
      parcs <- c('100', '400', '1000')
      
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
        png(paste0("./gta/connectomes/",mod, "_", dataset, "_", parc, ".png"), units = "px", width = 500, height = 500)
        image(Mplot, axes=FALSE, main=paste(mod, "mean, granularity:", parc), col=Col)
        dev.off()
        
        # ------------------------------------------------------------------------------- #
        # EDGES: Subject to group correlations
        # Correlation of the edges | upper triangle only
        Cor <- c()
        for (i in 1:N) {
          Cor <- c(cor(c(mtxs[,,i][upper.tri(mtxs[,,i])]), c(Mavg[upper.tri(Mavg)])), Cor)
        }
        
        # ------------------------------------------------------------------------------- #
        # Graph theoretical analysis: Subject to group correlations
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
        
        # Computational Time
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
      }
      
    } else {
      print(paste('[INFO].... NOTHING found connectomes in', dataset,mod))
    }
    
  }
  
  # save the Modality by dataset
  write.csv(corr, file=paste0(dataset, "_correlations.csv"), row.names = FALSE)
  write.csv(means, file=paste0(dataset, "_node-means.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------- #
#### Boxplot of the subject-group correspondance ####
datasets <- c('MICS', 'CamCAN', 'EpiC', 'EpiC2', 'sudmex', 'MSC', 'audiopath')
for (dataset in datasets) {
  data.corr <- c()
  try(data.corr <- read.csv(paste0(dataset,"_correlations.csv")))
  
  if (!is.null(data.corr)) {
    print(paste('[INFO].... plotting', dataset))
    data.corr$granularity <- as.factor(data.corr$granularity)
    
    # read the gradients by granularity and parcelation
    grains <- c(100, 400, 1000)
    modalities <- c('GD', 'SC', 'FC', 'MPC')
    G1 <- c()
    for (grain in grains) {
      for (mod in modalities) { Gs <- c()
        File <- paste0("csv_rho/",dataset,"-",as.character(grain),"_", mod, "-rho.csv")
        if (file.exists(File)) {
          Gs <- read.csv(File, header = FALSE)[,1]
          size <- length(Gs)
          G1 <- rbind(G1, data.frame(granularity=rep(grain, size), modality=rep(mod, size), G1=Gs))
        }
      }
    }
    
    # Gradients: Create a unique id for each modality/granularity
    id <- paste0(G1$granularity, '-', G1$modality)
    idN <- c(); for (i in table(id)) { idN <- c(idN, 1:i)}
    G1$id <- paste0(id, '-',sprintf("%03d", idN))

    # GTA: Create a unique id for each modality/granularity
    data.corr$id <- paste0(data.corr$granularity, '-', data.corr$modality)
    data.corr <- data.corr[order(data.corr$id),]
    idN <- c(); for (i in table(data.corr$id)) { idN <- c(idN, 1:i)}
    data.corr$id <- paste0(data.corr$id, '-',sprintf("%03d", idN))
    
    # merge data
    data.corr <- merge(data.corr, G1[,3:4], by = 'id')
    data.corr$granularity <- as.factor(data.corr$granularity)
    
    for (mod in modalities) {
      # filter data by modality
      Data <- data.corr %>%
        filter(modality==mod)
      lim <- 0
      
      # plot and save the data
      svg(filename=paste0("./gta/boxplots/",dataset,'-',mod,"_SubjectGroup-rho.svg"), width=16 , height=3 , pointsize=300, bg='transparent')
      g=grid.arrange(
        # Gradient 1
        ggplot(Data, aes(x=granularity, y=G1)) + 
          geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') + 
          xlab("Granularity") +
          ylim(lim, 1) +
          ggtitle(paste0(dataset, '-', mod, ': ', 'Gradient 1')) +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_jitter(width = 0.15, size = 1, colour='gray10', alpha=0.6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")),
        
        # Edges
        ggplot(Data, aes(x=granularity, y=edges)) + 
          geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') + 
          xlab("Granularity") +
          ylim(lim, 1) +
          ggtitle('Edges') +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_jitter(width = 0.15, size = 1, colour='gray10', alpha=0.6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")),
        
        # strength
        ggplot(Data, aes(x=granularity, y=strength)) + 
          geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') + 
          xlab("Granularity") +
          ylim(lim, 1) +
          ggtitle('Strength') +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_jitter(width = 0.15, size = 1, colour='gray10', alpha=0.6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")),
        
        # Charasteristic path length
        ggplot(Data, aes(x=granularity, y=path)) + 
          geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') + 
          xlab("Granularity") +
          ylim(lim, 1) +
          ggtitle('Char path length') +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_jitter(width = 0.15, size = 1, colour='gray10', alpha=0.6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")),
        
        # Weighted Cluster coefficient
        ggplot(Data, aes(x=granularity, y=cluster)) + 
          geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') + 
          xlab("Granularity") +
          ylim(lim, 1) +
          ggtitle('Cluster coefficient') +
          theme(plot.title = element_text(hjust = 0.5)) +
          geom_jitter(width = 0.15, size = 1, colour='gray10', alpha=0.6) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")),
        ncol=5, nrow = 1
      )
      plot(g)
      dev.off()
    }

    # Long to wide format    
    data.corr$id <- NULL
    Nmin <- min(table(data.corr$granularity))
    data.corr.wide <- data.corr[data.corr$granularity==100,]
    data.corr.wide <- data.corr.wide[1:Nmin,]
    for (i in c(400, 1000)) { data.corr.wide <- cbind(data.corr.wide, data.corr[data.corr$granularity==i,3:7][1:Nmin,]) }
    colnames(data.corr.wide)[3:17] <- paste0(rep(colnames(data.corr.wide)[3:7],3), rep(c('.100','.400','.1000'), each=5))
    data.corr.wide$granularity <- NULL
    
    # Correlation plot of all the variables 
    data.corr.mtx <- aggregate(cbind(G1.100,G1.400,G1.1000,edges.100,edges.400,edges.1000,strength.100,strength.400,strength.1000,path.100,path.400,path.1000,cluster.100,cluster.400,cluster.1000)~modality, data = data.corr.wide, mean)
    rownames(data.corr.mtx) <- data.corr.mtx$modality
    data.corr.mtx <- data.corr.mtx[c('GD', 'SC', 'FC', 'MPC'), ]
    data.corr.mtx$modality <- NULL
    data.corr.mtx[is.na(data.corr.mtx)] <- 0
    data.corr.mtx <- as.matrix(data.corr.mtx)
    

    # Plot and save the correlation matrix
    png(paste0(dataset, "_correlation.png"), units = "px", width = 700, height = 350)
    corrplot(data.corr.mtx, is.corr = FALSE, method = 'color',addCoef.col = 'gray10', title = dataset, cl.lim = c(0,1),cl.ratio = 0.6, cl.pos = 'b', col = rep(brewer.pal(9,'Reds'),2), tl.col = 'black')
    dev.off()
    
  } else {
    print(paste('[INFO].... EMPTY csv table:', dataset))
  }
}

# ------------------------------------------------------------------------------- #
#### Surface plot of the mean GTA ####
# fsverage5 read surfaces
pial.lh <- read.fs.surface(filepath = '/Users/rcruces/git_here/micapipe/surfaces/fsaverage5/surf/lh.pial')
pial.rh <- read.fs.surface(filepath = '/Users/rcruces/git_here/micapipe/surfaces/fsaverage5/surf/rh.pial')

# My color function
mymagma <- colorRampPalette(c(rep("gray75",2), brewer.pal(9, 'YlOrRd')))

# Plot Nodal features on surfaces by dataset, modality and granularity
for (dataset in c('MICS', 'CamCAN', 'EpiC', 'EpiC2', 'audiopath', 'sudmex', 'MSC')) {
  
  # Read dataset
  data.mean <- read.csv(paste0(dataset, "_node-means.csv"))
  
  # for each each modality
  for (mod in c('GD', 'SC', 'FC', 'MPC')) {
    
    # for each feature
    for (feat in c('strength', 'path', 'cluster')) {
      # for each granularity
      for (grain in c(100, 400, 1000)){
        # fsaverage5: read annotation files
        # left
        lh_annot_file <- paste0("/Users/rcruces/git_here/micapipe/parcellations/lh.schaefer-",as.character(grain),"_mics.annot")
        lh_annot = read.fs.annot(lh_annot_file)
        # right
        rh_annot_file <- paste0("/Users/rcruces/git_here/micapipe/parcellations/rh.schaefer-",as.character(grain),"_mics.annot")
        rh_annot = read.fs.annot(rh_annot_file)
        
        # labels
        atlas_region_names.lh <- lh_annot$colortable$struct_names
        atlas_region_names.rh <- rh_annot$colortable$struct_names  
        
        # Region-based results
        annot_values <- data.mean %>% filter(granularity==grain, modality==mod)
        if (length(annot_values$roi)>0) {
          val <- annot_values[,c(feat)]
          N <- length(val)
          
          # Normalize the values
          annot_values <- norm01(val)
          
          # left 
          annot_values.lh <- c(NaN, annot_values[1:(N/2)])
          names(annot_values.lh) <- atlas_region_names.lh
          # right
          annot_values.rh <- c(NaN, annot_values[((N/2)+1):N])
          names(annot_values.rh) <- atlas_region_names.rh
          
          # Set the color limits (1st and 3rd quantile)
          lf= limit_fun(quantile(annot_values, prob=0.2), quantile(annot_values, prob=0.95))
          lf = limit_fun(0,1)
          # Values to vertex
          val.lh <- spread.values.over.annot(annot = lh_annot, region_value_list = annot_values.lh, value_for_unlisted_regions = NaN)
          val.rh <- spread.values.over.annot(annot = rh_annot, region_value_list = annot_values.rh, value_for_unlisted_regions = NaN)
          
          # Create the coloredmeshes
          cml <- coloredmesh.from.preloaded.data(pial.lh, morph_data = lf(val.lh$spread_data), hemi = 'lh', makecmap_options = list('colFn'=mymagma))
          cmr <- coloredmesh.from.preloaded.data(pial.rh, morph_data = lf(val.rh$spread_data), hemi = 'rh', makecmap_options = list('colFn'=mymagma))
          
          # cretate the color mesh
          colormesh <- brainviews(views = 't4', coloredmeshes=list('lh'=cml, 'rh'=cmr), rglactions = list('trans_fun'=limit_fun(0, 1), 'no_vis'=T))
          # save the figure
          vis.export.from.coloredmeshes(colormesh, grid_like = FALSE, view_angles = c('sd_lateral_lh', 'sd_lateral_rh'), colorbar_legend=paste0(feat," ", grain),
                                        img_only = TRUE, horizontal=TRUE, output_img = paste0("./gta/surfaces/",dataset,"_",mod,"_",feat,"_", grain,".png"))
          # clean env
          while (rgl.cur() > 0) { rgl.close() }; file.remove(list.files(path = getwd(), pattern = 'fsbrain'))
        }
      }
    }
  }
}
