## Calculating Genetic Diversity

This file takes one input VCF and one metadata csv with individual names and population names. The script subsets the VCF to specified levels of missing data, subsets individuals to only individuals you want to keep, then uses hierfstat to calculate Ho, Hs, Fis and Fst for each population/group. 

The script outputs CSVs of the mean and standard deviations for all statistics, summarized over populations of interest.
