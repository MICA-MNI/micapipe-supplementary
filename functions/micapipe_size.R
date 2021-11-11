# R script
#
# Table of micapipe output size
# Display   values  are  in 1024 bytes (1 kilobyte, K) and round to ceiling
# 1M = 1024K
# 1G = 1024000K
#
# R version 3.6.3
#
# Created by RRC on October 2021 (the second year of the pademic)

# ------------------------------------------------------ # 
# Libraries
library("gtsummary")
library("tidyverse")

# ------------------------------------------------------ # 
# Set working directory to repository
setwd("/Users/rcruces/git_here/micapipe-supplementary/")

# ------------------------------------------------------ # 
# Helper functions
csv_2_tbl <- function(csv) {
  # Load the data
  data <- read.csv(paste0("data/", csv), as.is = TRUE)
  
  # Split path into variables
  data$id <- sapply(strsplit(as.character(data$path), split='/', fixed=TRUE), function(x) (x[1]))
  data$dir <- sapply(strsplit(as.character(data$path), split='/', fixed=TRUE), function(x) (x[2]))
  data$subdir <- sapply(strsplit(as.character(data$path), split='/', fixed=TRUE), function(x) (x[3]))
  
  # Order data
  Ndata <- rbind(
    # proc_rsfmri
    data %>% 
      select(bits, id, dir, subdir)  %>%
      filter(dir=="func", is.na(subdir)) %>% 
      select(bits, dir, id),
    
    # proc_dwi
    data %>% 
      select(bits, id, dir, subdir)  %>%
      filter(dir=="dwi", is.na(subdir)) %>% 
      select(bits, dir, id),
    
    # xfm
    data %>% 
      select(bits, id, dir, subdir)  %>%
      filter(dir=="xfm", is.na(subdir)) %>% 
      select(bits, dir, id),
    
    # QC
    data %>% 
      select(bits, id, dir, subdir)  %>%
      filter(dir=="QC", is.na(subdir)) %>% 
      select(bits, dir, id),
    
    # proc_structural | proc_freesurfer
    data %>% 
      select(bits, id, dir, subdir)  %>%
      filter(dir=="anat", is.na(subdir) ) %>% 
      select(bits, dir, id)
  )
  
  # Kilobytes to Megabytes
  Ndata$Kb <- round(as.numeric(as.character(Ndata$bits))/1024,0)
  
  # Plot the table
  nom <- gsub("micapipe_size_", "", gsub(".csv", "", csv))
  
  # Report percentage as in a NEJM article: n/N (p)
  tbl <- Ndata %>% 
    select(Kb, dir)  %>%
    tbl_summary(missing = "no", by=dir,
                statistic = list(all_continuous() ~ "{mean} ± {sd}", 
                                 all_categorical() ~ "{n} ({p}%)")) %>%
    modify_header(label = nom) %>%
    bold_labels() %>%
    italicize_levels() %>% as_gt()
  print(tbl)
}

# ------------------------------------------------------ # 
# Print the directory size per database
csvs <- list.files("data", pattern = "micapipe_size")
for (i in 1:length(csvs)) {
  csv_2_tbl(csvs[i])
}

