###  Chinook panel: Individual and Pop demographic plotting
###    Exploratory analysis to visualize the locations and population numbers
###    
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
Chinook_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", sep="", header = TRUE, colClasses="factor")
dim(Chinook_genepop)
Chinook_genepop[1:15,1:15]

POP_INFO <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS.txt", header =TRUE)
head(POP_INFO)
dim(POP_INFO)