# Calculate genetic diversity statistics for a variety of missing data levels 

# If you haven't already, be sure to open RStudio from Terminal using open /Applications/RStudio.app

library(hierfstat)
library(vcfR)
library(tidyverse)
library(ggtree)


# Set working directory where you want the results and new VCFs to go
setwd("~/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/genetic-diversity/") 

# Import pop data with individuals you want to keep in the VCF
pop.data <- read_delim("/Users/AirAlex/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/low-missing-inds-only/metadata-86-inds.csv", delim = ",")

# Set main VCF location
main.vcf <- "/Users/AirAlex/Documents/NorthCarolina/TangledBank/ast/data-analysis-combined-runs/111inds-no-confiscated-inds/88clust-111inds.vcf"

name.for.run <- "86inds"

# Set missing data thresholds
prop.missing.data <- c( 0.75, 0.8, 0.85)

# Start the calculations
for (i in prop.missing.data){
  
  ### Make VCF with new missing data threshold
  print(paste("Making VCF with ", i, " proportion missing data", sep = ""))
  
  # Make list of individuals to keep.
  # NB: Change sample to column with individual names
  write.table(pop.data$ind, "inds-to-keep.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Make the VCF
  system(paste("vcftools --vcf ", main.vcf, " --max-missing ", i, " --max-alleles 2 --min-alleles 2 --keep inds-to-keep.txt --recode-INFO-all --recode --out ", name.for.run, "-", i, "missing", sep = ""))
  
  new.vcf <- grep(pattern = "*.recode.vcf", list.files(), value = TRUE)
  
  # Subset to 1 SNP per locus
  system(paste("python --version")) # Note the next script must be run in python 2.7
  
  system(paste("python ~/Documents/NorthCarolina/TangledBank/tbc-hub/3rad-analyses/subsetSNPs.py ", new.vcf, " ", name.for.run,"-",i,"missing.vcf", sep = ""))
  
  system(paste("rm ", new.vcf))
  
  # Find just made vcf with single SNP per locus
  single.snp.search.term <- paste("*",i, "missing.vcf", sep="")
  single.snp.vcf <- grep(pattern = single.snp.search.term, list.files(), value = TRUE)
  
  # Read VCF
  print("Loading VCF...")
  my.vcf <- vcfR::read.vcfR(single.snp.vcf)
  
  # Convert to genind
  print("Converting to genind")
  genind <- vcfR2genind(my.vcf)
  
  # Match pop names to individuals
  inds <- data.frame(ind = row.names(genind$tab))
  pop.names <- left_join(inds, pop.data, by = "ind")$pop
  
  # Convert to hierfstat df
  print("Converting to HierFstat format...")
  df <- genind2hierfstat(genind, pop = pop.names)

  
  # CALCULATE!
  print("Calculating pop gen stats...")
  stats <- basic.stats(df)

  hs <- Hs(df)
  ho <- Ho(df)
  fis <- 1-(ho/hs)
 
  
  # WC84 Fst
  # print(paste("Calculating Fst for ", i, " proportion missing data. Started at ", Sys.time(), sep = ""))
  # wc84 <- as.matrix(genet.dist(df, method = "WC84"))

  ### Dosage based FST
  # Load vcf directly for dosage based FST
  # Load VCF to hierstat
  print("Loading VCF for FSTs..")
  myvcf <- read.VCF(single.snp.vcf)
  
  print(paste("Calculating Fst for ", i, " proportion missing data. Started at ", Sys.time(), sep = ""))
  fs.dat <- fs.dosage(myvcf, pop = pop.names)
  
  pairwise.fst <- fs.dat$Fst2x2
  
  
  
  # Write to disk
  print("Writing to disk...")

  write_csv(x = data.frame("pop" = names(ho), "ho" = ho), file = paste(name.for.run,"-",i,"missing-ho.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame("pop" = names(hs), "hs" = hs), file = paste(name.for.run,"-",i,"missing-hs.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame("pop" = names(fis), "fis" = fis), file = paste(name.for.run,"-",i,"missing-fis.csv", sep=""), col_names = TRUE)
  write_csv(x = data.frame(pairwise.fst), file = paste(name.for.run,"-",i,"missing-dosage-basedFST.csv", sep=""), col_names = TRUE)
  
  print(paste("Done with ", i, " proportion missing data!", sep = ""))
}








## Calculate the means +/- SD for each stat, across all the missing value thresholds
ho <- lapply(list.files(pattern = "*missing-ho.csv"), function(x) read.delim(x,sep = ",")) %>%
  bind_cols() %>%
  select(c("pop...1", contains("ho")))
hs <- lapply(list.files(pattern = "*missing-hs.csv"),function(x) read.delim(x,sep = ",")) %>%
  bind_cols() %>%
  select(c("pop...1", contains("hs")))
fis <- lapply(list.files(pattern = "*missing-fis.csv"), function(x) read.delim(x,sep = ",")) %>%
  bind_cols() %>%
  select(c("pop...1", contains("fis")))

# Fis
fis.mean <- data.frame(pop = fis[,1], mean = apply(fis[,2:ncol(fis)], 1, function(x) mean(x, na.rm= TRUE)))
fis.sd <- data.frame(pop = fis[,1], sd = apply(fis[,2:ncol(fis)], 1, function(x) sd(x, na.rm= TRUE)))


#### NOTE: CHANGE HERE IF POP != POPULATION UNIT
pop.counts <- count(pop.data, pop) %>%
  mutate(n = n, pop = pop)

fis.df <- left_join(fis.mean, fis.sd, by = "pop") %>%
  left_join(., pop.counts, by = "pop")

fis.final.df <- data.frame(pop = fis.df$pop, n = fis.df$n, fis = paste(round(fis.df$mean, 4), "±", round(fis.df$sd,2), sep = ""))
# 
write.table(fis.final.df, "fis_mean-and-sd_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Ho
ho.mean <- data.frame(pop = ho[,1], mean = apply(ho[,2:ncol(ho)], 1, function(x) mean(x, na.rm= TRUE)))
ho.sd <- data.frame(pop = ho[,1], sd = apply(ho[,2:ncol(ho)], 1, function(x) sd(x, na.rm= TRUE)))

ho.df <- left_join(ho.mean, ho.sd, by = "pop") %>%
  left_join(., pop.counts, by = "pop")

ho.final.df <- data.frame(pop = ho.df$pop, n = ho.df$n, ho = paste(round(ho.df$mean, 4), "±", round(ho.df$sd,2), sep = ""))
# 
write.table(ho.final.df, "ho_mean-and-sd_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Hs
hs.mean <- data.frame(pop = hs[,1], mean = apply(hs[,2:ncol(hs)], 1, function(x) mean(x, na.rm= TRUE)))
hs.sd <- data.frame(pop = hs[,1], sd = apply(hs[,2:ncol(hs)], 1, function(x) sd(x, na.rm= TRUE)))

hs.df <- left_join(hs.mean, hs.sd, by = "pop") %>%
  left_join(., pop.counts, by = "pop")

hs.final.df <- data.frame(pop = hs.df$pop, n = hs.df$n, hs = paste(round(hs.df$mean, 4), "±", round(hs.df$sd,2), sep = ""))
#
write.table(hs.final.df, "hs_mean-and-sd_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


## FST

# # Import all the FST matrices as a list, then remove the column with names (the names are the same as FIS/HO/HS)
fst <- lapply(list.files(pattern = "*dosage-basedFST.csv"), function(x) as.matrix(read.delim(x,sep = ",")))

# Calculate means and SDs
fst.mean <- simplify2array(fst) %>%
  apply(1:2, mean)
fst.sd <- simplify2array(fst) %>%
  apply(1:2, sd)

write.table(fst.mean, "fst_means_final.csv", sep = ",")

# Make cleaner version
fst.table.clean <- matrix(paste(round(fst.mean, 4), "±", round(fst.sd,2), sep = ""),nrow = nrow(fst.mean), ncol = ncol(fst.mean))
write.table(fst.table.clean, "fst_mean-and-sd_clean_final.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE )

# Convert mean and SD FST to migrants per generatios
mean.migrants <- apply(fst.mean, 1:2, function(x) (1-x)/(4*x))
sd.migrants <- apply(fst.sd, 1:2, function(x) (1-x)/(4*x))

migrants.table <- matrix(paste(round(mean.migrants, 2), "+/-", round(sd.migrants,2), sep = ""), nrow = nrow(mean.migrants), ncol = ncol(mean.migrants))
write.table(migrants.table, "migrants_table_mean-and-sd.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE )



  







