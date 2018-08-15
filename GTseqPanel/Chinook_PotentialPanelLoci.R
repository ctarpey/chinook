### To refine the potential panel markers for a refined Chinook panel 
###   That were chosen using a ranked FST of one snp per tag loci, combining the 3 data sets (Method 2)
###
###   Carolyn Tarpey | August 2018
### ---------------------------------------

### Potential panel loc
### This code requires the output of the Chinook_MarkerSelectionUsingFST_3DataSetsCombined.R 
### which has the FST and other info about the potential panel loci

### Population Maps 
### This code also requires a population map, which each individual of the analysis in the rows, and their population named in columns. 

### Previous panel identity
### Loci that have been part of a previous panel will be identified 

library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(hierfstat)
library(adegenet)
library(argparse)
library(stringi)
library(reshape2)
library(RColorBrewer)
#install.packages("argparse")

#______________________________________________________________________________________________
################################IMPORT YOUR DATA FILES

#input the text files that have the results for the runs of FST from Genepop for each of the 4 population groups, within each data 3 types. 

#allCombined RAD new loci data set estimated FST
potentialPanelmarkers_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/potentialPanelmarkers_FST.txt", header=TRUE,  row.names=NULL,
                                 stringsAsFactors = FALSE, na.strings = "-" )
dim(potentialPanelmarkers_FST_file)
head(potentialPanelmarkers_FST_file)

# ####################### May not need these, but they are the FST files for all the data sets, for all the pop groups:
# #allCombined RAD new loci data set estimated FST
# allCombined_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADallCombined_nonduplicated_FST_results.txt", header=TRUE,  row.names=NULL,
#                                  stringsAsFactors = FALSE, na.strings = "-" )
# dim(allCombined_FST_file)
# allCombined_FST_file[1:5,1:5]
# 
# ## RADGTseq data set estimated FST
# RADGTseq_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADGTSeq_nonduplicated_FST_results.txt", header=TRUE,  row.names=NULL,
#                               stringsAsFactors = FALSE, na.strings = "-" )
# dim(RADGTseq_FST_file)
# RADGTseq_FST_file[1:5,1:5]
# 
# ## GTseqONLY data set estimated FST
# GTseqONLY_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/GTSeqONLY_nonduplicated_FST_results.txt", header=TRUE, row.names=NULL,
#                                stringsAsFactors = FALSE, na.strings = "-" )
# dim(GTseqONLY_FST_file)
# GTseqONLY_FST_file[1:5,1:5]
# 
# ## new RAD Loci, AllmatchRADGTseq data set estimated FST
# AllmatchRADGTseq_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADallCombined_nondup_FST_results_ALLmatchRADGTSeq.txt", header=TRUE, row.names=NULL,
#                                       stringsAsFactors = FALSE, na.strings = "-" )
# dim(AllmatchRADGTseq_FST_file)
# AllmatchRADGTseq_FST_file[1:5,1:2]

####<-----------------------------------------------------------------------START HERE, figure out how to do the below
#Remove overlapping loci, but keep the version with the highest ranking. 

#figure out which panel the pre-existing panel loci were on. 


#figure out how many total were on the previous panels
potentialPanel_oldLoci <- potentialPanelmarkers_FST[potentialPanelmarkers_FST$Identity != "NewRAD",]
dim(potentialPanel_oldLoci)
potentialPanel_oldLociunique <- unique(potentialPanel_oldLoci$Locus)
length(potentialPanel_oldLociunique)

