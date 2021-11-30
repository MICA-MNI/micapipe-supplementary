# R script
#
# Table micapipe processing times
#       &
# Figure Sankey diagram
#
# R version 3.6.3
#
# Created by RRC on September 2021 (the second year of the pademic)

# ------------------------------------------------------ # 
# Libraries
library("gtsummary")
library("tidyverse")
library("igraph")
library("visNetwork")
library("networkD3")

# ------------------------------------------------------ # 
# Set working directory to repository
setwd("/Users/rcruces/git_here/micapipe-supplementary/")

# Print the processing time per database
csvs <- list.files("data", pattern = "micapipe_processed")

for (i in 1:length(csvs)) {
  data <- read_csv(paste0("data/", csvs[i]))
  nom <- gsub("micapipe_processed_sub_", "", gsub(".csv", "", csvs[i]))
  # Report percentage as in a NEJM article: n/N (p)
  tbl <- data %>% 
    select(Module, Processing.time)  %>%
    filter(Processing.time>0.9) %>% 
    tbl_summary(missing = "no", by=Module,
                       statistic = list(all_continuous() ~ "{mean} ± {sd}", 
                                        all_categorical() ~ "{n} ({p}%)")) %>%
    modify_header(label = nom) %>%
    bold_labels() %>%
    italicize_levels() %>% as_gt()
  print(tbl)
}


# ------------------------------------------------------ # 
#### micapipe steps by module #### 
mica <- read.csv("data/micapipe-steps.csv")
proc.rsfmri <- mica %>% 
  subset(module=="proc_rsfmri") 
G.rsfmri <-  graph_from_data_frame(proc.rsfmri, directed = F)

par(mar=c(1,1,1,1))
plot(G.rsfmri, edge.arrow.size = 0.1, vertex.shape="circle",
     vertex.label.color="black",vertex.size=2,
     vertex.label.cex=0.5,
     edge.color="gray80",
     vertex.label.degree=pi/2,
     layout=layout_as_tree(G.rsfmri, root = "BIDS_func", flip.y = TRUE))

# List all the layoputs
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite", layouts)]

# Iterate over each layout graph
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(G.rsfmri)) 
  plot(G.rsfmri, edge.arrow.size = 0.2, vertex.shape="circle", 
       vertex.label.color="black",
       vertex.label.cex=0.7,
       edge.color="gray80",
       layout=l, main=layout)}

# Interactive simple network
simpleNetwork(proc.rsfmri, height="300px", width="300px",        
              Source = 1,                 # column number of source
              Target = 2,                 # column number of target
              linkDistance = 10,          # distance between node. Increase this value to have more space between nodes
              charge = -100,                # numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
              fontSize = 14,               # size of the node names
              fontFamily = "serif",       # font og node names
              linkColour = "#666",        # colour of edges, MUST be a common colour for the whole graph
              nodeColour = "#69b3a2",     # colour of nodes, MUST be a common colour for the whole graph
              opacity = 0.9,              # opacity of nodes. 0=transparent. 1=no transparency
              zoom = T                    # Can you zoom on the figure?
)

# proc_structural
proc.struc <- mica %>%  subset(module=="proc_structural")
proc.struc.G <-  graph_from_data_frame(proc.struc, directed = T)
plot(proc.struc.G, edge.arrow.size = 0.1, vertex.shape="none",
     vertex.label.color="black",vertex.size=2,
     vertex.label.cex=0.75,
     edge.color="gray80",
     vertex.label.degree=pi/2,
     layout=layout_as_tree(proc.struc.G, root = "t1w", flip.y = TRUE))

simpleNetwork(proc.struc, height="300px", width="300px",        
              Source = 1,                 # column number of source
              Target = 2,                 # column number of target
              linkDistance = 10,          # distance between node. Increase this value to have more space between nodes
              charge = -200,                # numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
              fontSize = 14,               # size of the node names
              fontFamily = "serif",       # font og node names
              linkColour = "#666",        # colour of edges, MUST be a common colour for the whole graph
              nodeColour = "#69b3a2",     # colour of nodes, MUST be a common colour for the whole graph
              opacity = 0.9,              # opacity of nodes. 0=transparent. 1=no transparency
              zoom = T                    # Can you zoom on the figure?
)

# --------------------------------------------------- #
#### proc_structural Sankey diagram ####
proc.struc$from <- as.character(proc.struc$from)
proc.struc$to <- as.character(proc.struc$to)
Nodes <- data.frame(name=unique(c(proc.struc$from, proc.struc$to)))
Nodes$color <- rep("#A6A6A6", dim(Nodes)[1])
Nodes$group <- as.factor(Nodes$color)
Nodes$name <- as.character(Nodes$name)

Links <- proc.struc
colnames(Links)[1:2] <- c("source", "target")

for (i in 0:(dim(Nodes)[1]-1)) { 
  Links$source <- ifelse(Links$source==Nodes[i+1,1], i, Links$source)
  Links$target <- ifelse(Links$target==Nodes[i+1,1], i, Links$target) 
}
Links$group <- as.factor(ifelse(Links$value>0,1,0))
Links$source <- as.integer(Links$source)
Links$target <- as.integer(Links$target)
Links$value <-ifelse(Links$value>0,0.5,0)

# Links$value <-as.character(Links$command)
Links <- Links %>%  select(source, target, value, group)

# Sankey diagram of proc-structural
sankeyNetwork(Links = Links, Nodes = Nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", NodeGroup = "group", LinkGroup = "group",
              units = "", fontSize = 12, nodeWidth = 20, fontFamily = "Courier")

# --------------------------------------------------- #
#### proc_rsfmri Sankey diagram ####
proc.rsfmri$from <- as.character(proc.rsfmri$from)
proc.rsfmri$to <- as.character(proc.rsfmri$to)
Nodes <- data.frame(name=unique(c(proc.rsfmri$from, proc.rsfmri$to)))
Nodes$color <- rep("#104E8B", dim(Nodes)[1])
Nodes$group <- as.factor(Nodes$color)
Nodes$name <- as.character(Nodes$name)

Links <- proc.rsfmri
colnames(Links)[1:2] <- c("source", "target")

for (i in 0:(dim(Nodes)[1]-1)) { 
  Links$source <- ifelse(Links$source==Nodes[i+1,1], i, Links$source)
  Links$target <- ifelse(Links$target==Nodes[i+1,1], i, Links$target) 
}
Links$group <- as.factor(Links$proc)
Links$source <- as.integer(Links$source)
Links$target <- as.integer(Links$target)
Links$value <-ifelse(is.na(Links$value),0.5,0)
# Links$value <-as.character(Links$command)
Links <- Links %>%  select(source, target, value, group)

sankeyNetwork(Links = Links, Nodes = Nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", NodeGroup = "group", LinkGroup = "group",
              units = "", fontSize = 10, nodeWidth = 10, fontFamily = "Courier")

# --------------------------------------------------- #
#### proc_dwi Sankey diagram ####
proc.dwi <- mica %>% 
  subset(module=="proc_dwi") 

proc.dwi$from <- as.character(proc.dwi$from)
proc.dwi$to <- as.character(proc.dwi$to)
Nodes <- data.frame(name=unique(c(proc.dwi$from, proc.dwi$to)))
Nodes$color <- rep("#104E8B", dim(Nodes)[1])
Nodes$group <- as.factor(Nodes$color)
Nodes$name <- as.character(Nodes$name)

Links <- proc.dwi
colnames(Links)[1:2] <- c("source", "target")

for (i in 0:(dim(Nodes)[1]-1)) { 
  Links$source <- ifelse(Links$source==Nodes[i+1,1], i, Links$source)
  Links$target <- ifelse(Links$target==Nodes[i+1,1], i, Links$target) 
}
Links$group <- as.factor(Links$proc)
Links$source <- as.integer(Links$source)
Links$target <- as.integer(Links$target)
Links$value <-ifelse(is.na(Links$value),0.5,0)

# Links$value <-as.character(Links$command)
Links <- Links %>%  select(source, target, value, group)

sankeyNetwork(Links = Links, Nodes = Nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name", NodeGroup = "group", LinkGroup = "group",
              units = "", fontSize = 10, nodeWidth = 10, fontFamily = "Courier", nodePadding = 2)


