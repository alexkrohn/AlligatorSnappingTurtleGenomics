# Calculate genetic diversity statistics using hierfstat and ade4

# If you haven't already, be sure to open RStudio from Terminal using open /Applications/RStudio.app


library(hierfstat)
library(adegenet)
library(dplyr)
library(readr)
library(ggtree)
library(stringr)
library(progress)

## Convert VCFs to PLINK at a variety of missing datalevels
setwd("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/gen-diversity/")


missing.data.levels <- c(0.75, 0.8, 0.85, 0.9)
project.title <- "104inds"
vcf.location <- "~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/104inds-88clust.vcf"

convert.the.vcfs <- function(missing.data.levels, project.title, vcf.location){
  # Set progress bar
  for (i in missing.data.levels){
    system(paste("vcftools --vcf ", vcf.location, " --recode-INFO-all --recode --max-missing ", i, " --out ", project.title, str_remove(i, "\\."),"maxmissing", sep=""))
  }
  for (i in list.files(pattern = "*maxmissing.recode.vcf")){
    system(paste("plink -vcf ", i, " -recode --allow-extra-chr ", "-out ", str_remove(i, ".recode.vcf"), sep = ""))
  }
  for (i in list.files(pattern = "*.map")){
    system(paste("cut -f2 ", i, " | cut -f1 -d ':' > mycontigs.txt && cut -f2,3,4", i, " > mymap.txt && paste mycontigs.txt mymap.txt > ", i, " && rm mycontigs.txt mymap.txt", sep = ""))
  }
  for (i in str_remove(list.files(pattern = "*.map"), ".map")){
    system(paste("plink --chr-set 95 --allow-extra-chr --make-bed --recode A --out ", i, " --file ", i, sep=""))
  }
}

# Load the datasets, and calculate Hs, Ho, FIS and WC FST for each population. Output those calculations in a csv for each pop.
pb <- progress_bar$new(
  format = " Calculating pop gen stats [:bar] :percent :elapsedfull",
  total = 50, width= 60)

for (i in 1:length(list.files(pattern = "*.raw"))){
  # Load
  astgenlight <- read.PLINK(file = list.files(pattern = "*.raw")[i], map = list.files(pattern = "*.map")[i])
  
  pb$tick(0)
  # Convert to genind
  astgenind <- dartR::gl2gi(astgenlight, probar=TRUE,verbose = 5)
  
  pb$tick()
  # Convert to hierfstat df
  df <- genind2hierfstat(astgenind)
  
  pb$tick()
  # Remove astgenind and asgenlight to save memory
  astgenind <- NULL
  astgenlight <- NULL
  
  pb$tick()
  # Import pop data
  pop.data <- read_delim("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/104inds-just-inds-and-pops.csv", delim = ",")
  
  # Remove other columns
  pop.data <- dplyr::select(pop.data, ind, pop)
  
  # Change pop column to ind, and remove old ind name
  df <- rename(df, ind = pop)
  
  # bind df and pop data 
  df <- left_join(df, pop.data, by = "ind")
  
  final.df <- data.frame(pop = df$pop, dplyr::select(df, -pop, -ind))
  
  pb$tick()
  
  # CALCULATE!
  stats <- basic.stats(final.df)
  
  pb$tick()
  # Ho
  ho <- apply(stats$Ho,2,function(x){mean(x,na.rm=TRUE)})
  
  pb$tick()
  # Hs
  
  hs <- apply(stats$Hs,2,function(x){mean(x,na.rm=TRUE)})
  pb$tick()
  # Fis
  fis <- apply(stats$Fis,2,function(x){mean(x,na.rm=TRUE)})
  pb$tick()
  # WC84 Fst
  wc84 <- genet.dist(final.df, method = "WC84") 
  pb$tick()
  # file name
  file.name <- str_remove(list.files(pattern = "*.raw"), ".raw")[i]
  
  # Write to disk
  write_csv(x = data.frame("pops" = names(ho), "ho" = ho), file = paste(file.name, "-ho.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame("pops" = names(hs), "hs" = hs), file = paste(file.name, "-hs.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame("pops" = names(fis), "fis" = fis), file = paste(file.name, "-fis.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame("pops" = names(wc84), data.frame(as.matrix(wc84))), file = paste(file.name, "-fst-wc84.csv", sep=""), col_names = TRUE)
  pb$tick()
}
