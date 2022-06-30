## Attempt to quantify IBD!
setwd("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/ibd/")

# Plot pairwise genetic and geographic distance
library(RgoogleMaps)
library(raster)
library(gdata)
library(tidyverse)

# Use mean pairwise Fsts
fsts <- read.table("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/genetic-diversity/fst_means_final.csv", sep = ",", header = TRUE)

# Remove negative values (from pops with n=1 or 2)
fsts[fsts <= 0] <- NA
fsts[fsts == 0] <- NA

# Load lat, longs and population identifications. Note the lat/long info is not available on GitHub due to poaching concerns.
metadata <- read.table("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/metadata-86-inds.csv", sep = ",", header = TRUE, fill = TRUE)

# Take mean lat and longs for each population
mean.pop.latlongs <- metadata %>%
  group_by(pop) %>%
  summarise(lat = mean(Latitude), long = mean(Longitude))

# Set coordinate system
coordinates(mean.pop.latlongs)<-c("long","lat")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
proj4string(mean.pop.latlongs) <- crs.geo

# Create geographic distance matrix in meters
geodist<-pointDistance(mean.pop.latlongs,longlat=TRUE)
geodist[which(geodist==0)] <- NA


# Create a matrix that describes how to color the dots: one color for all Suwannee comparisons, one for all Appalachicolae comparisons and a third for all temminckii comparisons
# Label all Suwannee Comparisons
suwannee.position <- mean.pop.latlongs$pop == "Suwannee"
suwannee.matrix <- outer(suwannee.position, suwannee.position, function(a, b) as.integer(a + b == 1L))
suwanee.matrix[suwanee.matrix == 1] <- "dodgerblue"

# Label all Apalachicolae comparisons
ap.positions <- mean.pop.latlongs$pop == "Apalachicola" | mean.pop.latlongs$pop == "Econfina" | mean.pop.latlongs$pop == "Chipola" | mean.pop.latlongs$pop == "Ochlockonee"
ap.matrix <- outer(ap.positions, ap.positions, function(a, b) as.integer(a + b == 1L))
ap.matrix[ap.matrix == 1] <- "tomato2"

# Plot each comparison 
par(mar=c(5.1,4.7,4.2,2.1))
# Apalachicolae first
plot(geodist[ap.matrix == 1]/1000, as.matrix(fsts[ap.matrix == 1]), ylab = expression('Pairwise Genetic Distance (F'[ST]*')'), xlab = "Pairwise Geographic Distance (km)", pch = 19, col = "tomato2", cex=1.5, cex.lab = 1.3, xlim = c(0, max(geodist/1000, na.rm=TRUE)), ylim = c(0, max(fsts, na.rm = TRUE)))
points(x = geodist[suwannee.matrix == 1]/1000,y = as.matrix(fsts[suwannee.matrix == 1]), pch = 19, col = "dodgerblue", cex=1.5, cex.lab = 1.3)
points(x = geodist[suwannee.matrix == 0 & ap.matrix == 0]/1000,y = as.matrix(fsts[suwannee.matrix == 0 & ap.matrix == 0]), pch = 19, col = "black", cex=1.5, cex.lab = 1.3)
legend(x=954,y = 0.12,col=c("dodgerblue2","tomato2","black"),pch=19,legend=c("Eastern Assemblage", "Central Assemblage", "Other Comparisons"), bg="transparent", xjust = 0)




# Use Mantel Test to evaluate whether river drainage difference or geographic distance better explains divergence

# Add lower trinagle values to upper triangle of geodist
geodist[upper.tri(geodist)] <- geodist[lower.tri(geodist)]

# Remove rows/cols with all NAs (Blackwater, Chipola, Jourdan)
fsts.no.nas <- fsts[rowSums(is.na(fsts)) != ncol(fsts),rowSums(is.na(fsts)) != ncol(fsts)]
geodist.no.nas <- geodist[rowSums(is.na(fsts)) != ncol(fsts),rowSums(is.na(fsts)) != ncol(fsts)]

# Test for the effect of IBD
vegan::mantel(geodist.no.nas,fsts.no.nas, permutations=999)

# Call:
#   vegan::mantel(xdis = geodist.no.nas, ydis = fsts.no.nas, permutations = 999) 
# 
# Mantel statistic r: 0.5221 
# Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
#   0.123 0.158 0.192 0.235 
# Permutation: free
# Number of permutations: 999

