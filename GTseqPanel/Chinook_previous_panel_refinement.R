###  Chinook panel: filtering existing panel loci for a refined GT-seq panel
###    Using the combined Rad/taq/Amp data file Garrett made
###    these are 847 loci that were included in the previous panels, re-examined
### Carolyn Tarpey | June 2018 
### ---------------------------------------


library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(reshape2)
library(plotly)
library(gridExtra)
library(scales) 
library(grid)


#load the genepop file with all the individuals and loci:  
Chin_847_geno_loci <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", sep="", header = TRUE, colClasses="factor")
dim(Chin_847_geno_loci)
Chin_847_geno_loci[1:15,1:15]

POP_INFO <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS.txt", header =TRUE)
head(POP_INFO)
dim(POP_INFO)
