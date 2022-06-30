## Attempt to quantify IBD!
setwd("AlligatorSnappingTurtleGenomics/ibd/")

# Plot pairwise genetic and geographic distance
library(RgoogleMaps)
library(raster)
library(gdata)
library(dplyr)

# Use mean pairwise Fsts (output from calc-genetic-diversity-for-a-variety-of-missing-data)
fsts <- read.table("AlligatorSnappingTurtleGenomics/genetic-diversity/mean-fsts.csv", sep = ",", header = TRUE)

# Remove negative values (from pops with n=1 or 2)
fsts[fsts <= 0] <- NA
fsts[fsts == 0] <- NA

# Load lat longs. Note, this metadata is not available due to poaching concerns. It is simply a csv with latitude, longitude and individual names for each sample.
lats <- read.table("metadata/86-inds-with-latlong.csv", sep = ",", header = TRUE, fill = TRUE)

# Join with the pop dataframe to assign individuals to populations, then take the mean of each populations lat and long values
pop.info <- read.table("../metadata/metadata-86-inds.csv", header=TRUE, sep = ",")

mean.pop.latlongs <- lats %>%
  left_join(pop.info, by = "ind") %>%
  group_by(pop.x) %>%
  summarise(lat = mean(Latitude), long = mean(Longitude))

# Set coordinate system
coordinates(mean.pop.latlongs)<-c("long","lat")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(mean.pop.latlongs) <- crs.geo

# Create geographic distance matrix
geodist<-pointDistance(mean.pop.latlongs,longlat=TRUE)
geodist[which(geodist==0)] <- NA

# # Color by whether the comparison involves an individual from Suwannee (1), Appalachicolae (2), or temenicki (3)
color.matrix <- as.matrix(read.table("color-by-population-comparisons.csv", sep = ",", header = TRUE, row.names = 1))

color.matrix[color.matrix == 1] <- "dodgerblue"
color.matrix[color.matrix == 2] <- "tomato2"
color.matrix[color.matrix == 3] <- "black"

# Make sure they're only taking the lower triangles
color.matrix[lower.tri(color.matrix, diag = FALSE) == FALSE] <- NA
geodist[lower.tri(geodist, diag = FALSE) == FALSE] <- NA
fsts[lower.tri(fsts, diag= FALSE) == FALSE] <- NA

# Plot
par(mar=c(5.1,4.7,4.2,2.1))
plot(geodist/1000, as.matrix(fsts), ylab = expression('Pairwise Genetic Distance (F'[ST]*')'), xlab = "Pairwise Geographic Distance (km)", pch = 19, col = color.matrix, cex=1.5, cex.lab = 1.3)
legend(x=954,y = 0.12,col=c("dodgerblue2","tomato2","black"),pch=19,legend=c("Eastern Assemblage", "Central Assemblage", "Other Comparisons"), bg="transparent")


# Use Mantel Test to evaluate whether river drainage difference or geographic distance better explains divergence

# First remove rows/cols with NAs (Wetappo/Jourdan)
geodist.no.nas <- geodist[c(1:10,12:17,19:24,26),c(1:10,12:17,19:24,26)]
fsts.no.nas <- fsts[c(1:10,12:17,19:24,26),c(1:10,12:17,19:24,26)]

as.matrix(fsts)

# Test for the effect of IBD
vegan::mantel(geodist.no.nas,fsts.no.nas, permutations=999)
# Call:
#   vegan::mantel(xdis = geodist.no.nas, ydis = fsts.no.nas, permutations = 999) 
# 
# Mantel statistic r: 0.5133
# Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#   0.122 0.156 0.193 0.219 
# Permutation: free
# Number of permutations: 999
