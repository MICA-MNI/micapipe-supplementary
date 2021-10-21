# R script
#
# Suplementary figure micapipe 
# Subject-Group mean gradient Correlation coefficients
# R version 3.6.3
#
# Created by RRC on September 2021 (the second year of the pademic)

# libraries
library(tidyverse)
library(svglite)

# Set the working dir
setwd("/Users/rcruces/git_here/micapipe-supplementary/data/csv")

# Iterative variables
mod <- c("_FC", "_SC", "_GD", "_MPC")
lim <- c(0, 0.5, 0.8, 0)

# For each Modality / file
for (k in 1:4) {
  cnn <- mod[k]
  # List the files
  files <- list.files(path = ".", pattern = cnn)
  for (i in 1:length(files)) {
    # Load and sort the dataset rhos
    dataset <- strsplit(files[i], '_')[[1]][1]
    Data <- read.csv(files[i], header = FALSE, sep = ',')
    colnames(Data) <- paste0("G", sprintf('%0.2d', 1:10))
    N <- dim(Data)[1]
    Data <- gather(Data, "gradients", "rho",1:10)
    
    # print info
    print(paste0(dataset,cnn,".... ", files[i]))
    
    # plot and save the data
    svg(filename=paste0("../svg/",dataset,cnn,"_group-rho.svg"), width=3.75 , height=2.75 , pointsize=300, bg='transparent')
    g=ggplot(Data, aes(x = gradients, y = rho)) +
        geom_jitter(width = .1, size = 1, colour='gray10', alpha=0.2) +
        ylim(lim[k], 1) +
        ggtitle(paste0(dataset,cnn, ", N=",N) )+
        geom_boxplot(fill="gray75", size = 0.35, outlier.shape = NA, alpha=1, colour='gray10') +
        theme(plot.title = element_text(hjust = 0.5))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    plot(g)
    dev.off()
  }
}
