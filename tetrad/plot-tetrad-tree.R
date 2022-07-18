# Try plotting tetrad and tetrad cluster tree

library(tidyverse)
library(ggtree)
library(treeio)

# Read tree
tetrad.tree <- read.newick("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/tree/tetrad/ast-tetrad-86inds.tree.cons", node.label = "support")

# Root at Suwannee
tetrad.tree.rooted <- ape::root(phy = tetrad.tree, outgroup = "MATE423", resolve.root = FALSE ) 

# Convert to tibble, then add correct names 
new.names <- read.table("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/metadata-86-inds.csv", header=TRUE, sep = ",") %>%
  mutate(new.names = paste(ind,"-", pop, sep = "")) %>%
  rename(label = ind) %>%
  dplyr::select(label, new.names)

rooted.tree.with.names <- as_tibble(tetrad.tree.rooted) %>%
  left_join(., new.names, by = "label")


## Plot the tree, but showing node labels on the nodes, so we know where to draw clades
ggtree(as.treedata(rooted.tree.with.names), branch.length = "none")+
  geom_tiplab(aes(label=new.names), size=2.75, hjust=-0.5)+
  geom_text(aes(label=node), hjust=-0.02)+
  xlim(0,30)

# Set colors based on the nodes/tips
# First, set all colors to black
rooted.tree.with.names$color <- rep("black", nrow(rooted.tree.with.names))

# Set Western Assemblage West of MS River as purple
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 40 & rooted.tree.with.names$node <= 86 | rooted.tree.with.names$node >= 125 & rooted.tree.with.names$node <= 170)] <- "mediumpurple2"

# Set FL/AL/MS Western Assemblage color
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 143 & rooted.tree.with.names$node <= 145 | rooted.tree.with.names$node >= 58 & rooted.tree.with.names$node <= 60 | rooted.tree.with.names$node >= 1 & rooted.tree.with.names$node <= 25 | rooted.tree.with.names$node >= 88 & rooted.tree.with.names$node <= 112 )] <- "sienna1"

# Set Apalachicolae color
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 28 & rooted.tree.with.names$node <= 39 | rooted.tree.with.names$node >= 113 & rooted.tree.with.names$node <= 124 )] <- "dodgerblue3"


# Plot tree with color
tiff("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/tree/tetrad/tetrad-tree-color.tiff", width = 5.2, height = 8.1, res = 350, units = "in")
ggtree(as.treedata(rooted.tree.with.names), branch.length = "none", aes(color = I(color)))+
  geom_tiplab(aes(label=new.names), size = 2.75)+
  geom_cladelab(node = 106, label = "Western\nAssemblage", offset = 11, extend = 0.5)+
  geom_cladelab(node = 114, label = "Central Assemblage", offset = 12, extend = 0.5)+
  geom_cladelab(node = 26, label = "Eastern Assemblage", offset = 10, offset.text = 0.005, extend = 1.4)+
  geom_nodepoint(aes(label = support, 
                     subset = support >= 85))+
  xlim(0,39)
dev.off()

# Plot in black and white
tiff("~/Documents/NorthCarolina/TangledBank/ast/writing-and-results/sena-major-revisions/figs/Fig6-msc-tree-black-and-white.tiff", width = 5.2, height = 8.1, res = 650, units = "in")
ggtree(as.treedata(rooted.tree.with.names), branch.length = "none")+
  geom_tiplab(aes(label=new.names), size = 2.75)+
  geom_cladelab(node = 106, label = "Western\nAssemblage", offset = 50, offset.text = 1)+
  geom_cladelab(node = 96, label = "Western Assemblage\nE of MS River", offset = 21, extend = 11, offset.text = 1)+
  geom_cladelab(node = 144, label = "Western Assemblage\nE of MS River", offset = 20, extend = 0, offset.text = 1)+
  geom_cladelab(node = 146, label = "Western Assemblage\nfrom MS River and TX", offset = 20, extend = 0, offset.text = 1)+
  geom_cladelab(node = 126, label = "Western Assemblage\nfrom MS River and TX", offset = 20, extend = 0, offset.text = 1)+
  geom_cladelab(node = 114, label = "Central Assemblage", offset = 22, extend = 0.5, offset.text = 1)+
  geom_cladelab(node = 26, label = "Eastern Assemblage", offset = 20, offset.text = 0.005, extend = 1.4, offset.text = 1)+
  geom_nodepoint(aes(label = support, 
                     subset = support >= 85))+
  xlim(0,80)
dev.off()
