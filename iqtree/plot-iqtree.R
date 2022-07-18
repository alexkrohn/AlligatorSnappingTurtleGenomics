## Plot trees from iqTree

## 86 individuals
library(tidyverse)
library(ggtree)
library(ape)
library(treeio)

# Consensus tree from iqtree
iqtree <- read.newick("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/tree/iqtree/iqtree-8maxmissing-biallelic-86inds.min1.phy.contree", node.label = "support")

# Root around Suwaneensis
iqtree.rooted <- ape::root(phy = iqtree, outgroup = "MATE423", resolve.root = FALSE ) 


# Convert to tibble, then add correct names 
new.names <- read.table("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/metadata-86-inds.csv", header=TRUE, sep = ",") %>%
  mutate(new.names = paste(ind,"-", pop, sep = "")) %>%
  rename(label = ind) %>%
  dplyr::select(label, new.names)

rooted.tree.with.names <- as_tibble(iqtree.rooted) %>%
  left_join(., new.names, by = "label")



# Convert back to tree
tree.with.names <- as.treedata(rooted.tree.with.names)




# Plot tree, with nodes bootstraps > 85 labelled
ggtree(tree.with.names)+
  geom_tiplab(aes(label=new.names), size=2.75)+
  geom_nodepoint(aes(label = support, 
                     subset = support >= 85))+
  xlim(0,1)+
  theme_tree2()


# Color the tree based on the PCA colors

# Create tree with no branch lengths, and with nodes labeled with their node number
ggtree(tree.with.names, branch.length = "none")+
  geom_tiplab(aes(label=new.names), size=2.75, hjust=-0.5)+
  geom_text(aes(label=node), hjust=-0.02)+
  xlim(0,30)


# Set colors based on the nodes/tips
# First, set all colors to black
rooted.tree.with.names$color <- rep("black", nrow(rooted.tree.with.names))

# Set Western Assemblage West of MS River as purple
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 89 & rooted.tree.with.names$node <= 114 | rooted.tree.with.names$node >= 1 & rooted.tree.with.names$node <= 44 | rooted.tree.with.names$node >= 115 & rooted.tree.with.names$node <= 131 & rooted.tree.with.names$node != 132 )] <- "mediumpurple2"

# Set FL/AL/MS Western Assemblage color
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 45 & rooted.tree.with.names$node <= 48 | rooted.tree.with.names$node >= 151 & rooted.tree.with.names$node <= 170 | rooted.tree.with.names$node >= 63 & rooted.tree.with.names$node <= 86 | rooted.tree.with.names$node == 88 | rooted.tree.with.names$node >= 132 & rooted.tree.with.names$node <= 138 )] <- "sienna1"

# Set Apalachicolae color
rooted.tree.with.names$color[which(rooted.tree.with.names$node >= 140 & rooted.tree.with.names$node <= 150 | rooted.tree.with.names$node >= 49 & rooted.tree.with.names$node <= 60 )] <- "dodgerblue3"


# Plot tree with color
tiff("~/Documents/NorthCarolina/TangledBank/ast/writing-and-results/sena-major-revisions/figs/Fig5-ml-tree.tiff", width = 5.2, height = 8.1, unit = "in", res = 350)
ggtree(as.treedata(rooted.tree.with.names), aes(color = I(color)))+
  geom_tiplab(aes(label=new.names), size = 2.75)+
  geom_nodepoint(aes(label = support, 
                     subset = support >= 85))+
  geom_cladelab(node = 138, label = "Western\nAssemblage", offset = 0.32)+
  geom_cladelab(node = 140, label = "Central Assemblage", offset = 0.35)+
  geom_cladelab(node = 61, label = "Eastern Assemblage", offset = 0.32, offset.text = 0.005, extend = 1.2)+
  xlim(0,1.33)
dev.off()


## Plot tree with black and white labels
tiff("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/tree/iqtree/mltree-scaled-for-paper-black-and-white.tiff", width = 5.2, height = 8.1, unit = "in", res = 650)
ggtree(tree.with.names)+
  geom_tiplab(aes(label=new.names), size = 2.75)+
  geom_nodepoint(aes(label = support, 
                     subset = support >= 85))+
  geom_cladelab(node = 138, label = "Western\nAssemblage", offset = 1.9)+
  geom_cladelab(node = 89, label = "Western\nAssemblage\nfrom MS River and TX", offset = 0.8, extend = 0.1)+
  geom_cladelab(node = 160, label = "Western Assemblage\nE of MS River", offset = 0.85, extend = 11)+
  geom_cladelab(node = 140, label = "Central Assemblage", offset = 0.85, extend = 0.3)+
  geom_cladelab(node = 61, label = "Eastern Assemblage", offset = 0.8, offset.text = 0.005, extend = 1.25)+
  xlim(0,3.2)
dev.off()

