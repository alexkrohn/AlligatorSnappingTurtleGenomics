## Calculating Genetic Diversity

This file takes one input VCF and one metadata csv with individual names and population names. The script subsets the VCF to specified levels of missing data, subsets individuals to only individuals you want to keep, then uses hierfstat to calculate H<sub>O</sub>, H<sub>S</sub>, F<sub>IS</sub> and F<sub>ST</sub> for each population/group. 

The script outputs CSVs of the mean and standard deviations for all statistics, summarized over populations of interest.
