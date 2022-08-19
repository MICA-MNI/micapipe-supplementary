#--------------------------------------------------- 
# Check bvecs rotations with MRtrix/ants
# Made for micapipe review
# Raul R. Cruces August, 17, 2022
#--------------------------------------------------- 

# Libraries
library("plot3D")
library("dplyr")
library("scales")

# Load
bvecs <- read.csv("bvecs.csv")

# Add origin columns to bvecs
bvecs$x0 <- rep(0,dim(bvecs)[1])   

#--------------------------------------------------- 
#### ORIGINAL BVECS ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = "red4",
           lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
           main = "Original bvecs", bty ="g", ticktype = "detailed")
 # Add starting point of arrow
points3D(0, 0, 0, add = TRUE, col="red4", colkey = FALSE, pch = 19, cex = 1)

#--------------------------------------------------- 
####  BVECS after MRTRIX handle rotation 45 ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec.x45", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = "coral2",
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "Rotated bvecs", bty ="g", ticktype = "detailed")
# Add starting point of arrow
points3D(0, 0, 0, add = TRUE, col="coral2", colkey = FALSE, pch = 19, cex = 1)

#--------------------------------------------------- 
####  ORIGINAL BVECS and MRtrix3 rotated 45 ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec" | group == "bvec.x45", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = c(alpha("red4",0.5), alpha("coral2",0.5))[subgr$group],
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "Original & Rotated bvecs", bty ="g", ticktype = "detailed")

#--------------------------------------------------- 
#### MRtriX3 mrconvert to mif: ORIGINAL ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "orig", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = "cornflowerblue",
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "MRtrix3 dwgrad: original bvecs", bty ="g", ticktype = "detailed")
# Add starting point of arrow
points3D(0, 0, 0, add = TRUE, col="cornflowerblue", colkey = FALSE, pch = 19, cex = 1)
#--------------------------------------------------- 
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec" | group == "orig", !num == 0)
subgr$group <- droplevels(subgr$group)
# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = c(alpha("red4",0.3), "cornflowerblue")[subgr$group],
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "Original & MRtrix3 dwgrad bvecs", bty ="g", ticktype = "detailed")
text3D(subgr$x, subgr$y, subgr$z, as.character(subgr$num),
       col = c(alpha("red4",0.5), "cornflowerblue")[subgr$group], add=TRUE, colkey = FALSE)

#--------------------------------------------------- 
#### MRtriX3 mrconvert to mif: mrtransform rot45 ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "x45.mrtx", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = "dodgerblue3",
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "MRtrix3 dwgrad: 45°bvecs", bty ="g", ticktype = "detailed")
# Add starting point of arrow
points3D(0, 0, 0, add = TRUE, col="dodgerblue3", colkey = FALSE, pch = 19, cex = 1)
#--------------------------------------------------- 
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec" | group == "x45.mrtx", !num == 0)
subgr$group <- droplevels(subgr$group)
# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = c(alpha("red4",0.3), "dodgerblue3")[subgr$group],
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "Original & 45°bvecs", bty ="g", ticktype = "detailed")
text3D(subgr$x, subgr$y, subgr$z, as.character(subgr$num),
       col = c(alpha("red4",0.5), "dodgerblue3")[subgr$group], add=TRUE, colkey = FALSE)

#--------------------------------------------------- 
#### MRtriX3 mrconvert to mif: antsapplytransform rot45 ####
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "x45.ants", !num == 0)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = "royalblue4",
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "MRtrix3 dwgrad: 45°ants bvecs", bty ="g", ticktype = "detailed")
# Add starting point of arrow
points3D(0, 0, 0, add = TRUE, col="royalblue4", colkey = FALSE, pch = 19, cex = 1)

#--------------------------------------------------- 
subgr <- bvecs %>%
  select(num, x0, x, y, z, group) %>%
  filter(group == "bvec" | group == "x45.ants", !num == 0)
subgr$group <- droplevels(subgr$group)

# 3D arrows plot from origin
arrows3D(subgr$x0, subgr$x0, subgr$x0, subgr$x, subgr$y, subgr$z, col = c(alpha("red4",0.5), "royalblue4")[subgr$group],
         lwd = 2, d = 3, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         main = "Original & 45°ants bvecs", bty ="g", ticktype = "detailed")
text3D(subgr$x, subgr$y, subgr$z, as.character(subgr$num),
      col = c(alpha("red4",0.5), "royalblue4")[subgr$group], add=TRUE, colkey = FALSE)
