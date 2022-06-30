# Make PCA from ipyrad PCA output

# On 104 Alligator Snapping Turtle individuals, 4010 SNPs, one SNP per locus post filtering
library(tidyverse)
library(ggnewscale)

# Import ipyrad PCA output of Eigen vectors
pca.data <- read.table("pca-df-mincov5-no-pops-with-1ind.csv", sep=",", skip = 1)

# Rename first column as ind to match csv
pca.data$ind <- pca.data$V1

# Import pop data
pop.data <- read.table("../metadata/104inds-just-inds-and-pops.csv", sep = ",", header=TRUE)

# Join them by ind
fulldataset <- left_join(pca.data, pop.data, by= "ind")

# Separate colors by Central, Western and Eastern Lineages with ggnewscale
library(ggnewscale)

## Add column identifying each individual into an assemblage
fulldataset.with.assemblage <- fulldataset %>%
  mutate(assemblage = case_when(pop == "Apalachicola" | pop == "Ochlockonee" | pop == "Chipola" | pop == "Econfina" | pop == "Wetappo" ~ "Central Assemblage",
                               pop == "Suwannee" ~ "Eastern Assemblage",
                             V2 < -9.5 ~ "Western Assemblage from MS River Westward",
                             TRUE ~ "FL/AL/MS Western Assemblage")) 

## Split the dataset into a list, separated by assemblage
assemblage.list <- split(fulldataset.with.assemblage, fulldataset.with.assemblage$assemblage)

# Set colors for each assemblage
eastern.cols <- "black"

# 5 colors for Central Assemblage
central.cols <- c("blue4", "royalblue2","steelblue1", "lightblue3",  "skyblue")

# 9 colors for Western Assemablage West of MS River
west.of.ms.cols <- c("purple4", "darkmagenta","mediumpurple3", "purple1","magenta2","orchid2", "darkorchid1", "plum1", "maroon1")

# 11 colors for the FL/AL/MS Western Assemblage rivers
fl.al.ms.western.cols <- c("darkred", "red3", "brown2", "tomato3", "coral2", "darkgoldenrod2", "darkorange2", "salmon1", "lightcoral", "tomato",  "sandybrown")

## Make the actual plot
ggplot(mapping = aes(x = V2, y = V3)) +
  geom_point(data = assemblage.list$`Eastern Assemblage`, aes(col = pop), size = 4) +
  scale_colour_manual(values = eastern.cols, 
                      guide = guide_legend(order = 1, title = expression(italic("M. suwanniensis"))))+
  new_scale_color() +
  geom_point(data = assemblage.list$`Central Assemblage`, aes(col = pop), size = 4) +
  scale_colour_manual(values = central.cols,
                      guide = guide_legend(order = 2, title = "Central Assemblage"))+
  new_scale_color()+
  geom_point(data = assemblage.list$`FL/AL/MS Western Assemblage`, aes(col = pop), size = 4) +
  scale_colour_manual(values = fl.al.ms.western.cols,
                      guide = guide_legend(order = 3, title = "FL/AL/MS Western Assemblage"))+
  new_scale_color()+
  geom_point(data = assemblage.list$`Western Assemblage from MS River Westward`, aes(col = pop), size = 4) +
  scale_colour_manual(values = west.of.ms.cols,
                      guide = guide_legend(order = 3, title = "Western Assemblage West of MS River"))+
  xlab("PC2 (11.5% variance explained)")+
  ylab("PC3 (9.9% variance explained)")+
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        legend.box.just = "left")
  

filter(fulldataset, V2 > 3.5 & V2 < 10)
