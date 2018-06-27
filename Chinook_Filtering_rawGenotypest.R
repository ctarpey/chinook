###  Chinook panel: filtering raw RAD genotypes 
###    for new RAD loci to included in a refined panel
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


#load the raw genepop file, edited for R. Includes all individuals and loci:  
Chin_raw_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/batch_13.allCombined_NoCookInlet_R_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(Chin_raw_genepop)
Chin_raw_genepop[1:15,1:15]
rownames(Chin_raw_genepop)<-gsub(",","",rownames(Chin_raw_genepop))
colnames(Chin_raw_genepop)<-gsub("X","",colnames(Chin_raw_genepop))


#a population map that has all the individuals and the populations that they belong to. 
Chin_all_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/populationMap_NoCookInlet.txt", header =FALSE)
colnames(Chin_all_popINFO) <- c("SampleName", "Population", "VeryBroadRegion", "Location", "BroadRegion", "FineReportingGroup")
head(Chin_all_popINFO)
dim(Chin_all_popINFO)

