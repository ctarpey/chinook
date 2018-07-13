###  Chinook panel: creating Genepop files 
###    for new RAD loci to included in a refined panel
###    for RADtaqAMP data set
###    for GTseqOnly data
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

### Read in the data files we will need ########

####Write a list of the individuals that have passed all the filters so far:
filteredSample_names <- read.table("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/filteredSample_names.txt", header=FALSE)
head(filteredSample_names)

###RAD LOCI
#load the filtered RAD loci genepop file, edited in R. #this file already has a population map merged in it
chin_RAD_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/filteredGenos_filteredSamples_filtpara.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_RAD_genepop)

#a population map that has all the individuals and the populations that they belong to. 
chin_RAD_all_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/populationMap_NoCookInlet.txt", header =TRUE)
head(chin_RAD_all_popINFO)
dim(chin_RAD_all_popINFO)

###RADtaqman amplicon LOCI
#load the RADtaqman amplicon LOCI genepop file
RADtaqAmp_Rcombined <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", header =TRUE, row.names= NULL, colClasses="factor")
RADtaqAmp_Rcombined[1:5,1:5]
RADtaqAmp_Rcombined[,849]<-NULL #there is an extra tab on the end of each line that is making an empty column 
dim(RADtaqAmp_Rcombined)

#a population map that has all the individuals and the populations that they belong to. 
RADtaqAmp_popINFO <- read.table("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/PopulationMap_847RAD_taq_genepop.txt", header=TRUE)
head(RADtaqAmp_popINFO)
dim(RADtaqAmp_popINFO)

###GTseq only loci
#load the RADtaqman amplicon LOCI genepop file
GTseq_ONLY_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/combinedPrimerProbe_amplicon_polyResults_genepop.txt", header =TRUE, colClasses="factor")
GTseq_ONLY_genepop[1:5,1:5]
dim(GTseq_ONLY_genepop)

#a population map that has all the individuals and the populations that they belong to. 
GTseq_ONLY_popINFO <- read.table("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/PopulationMap_GTseqDataONLY.txt", header=TRUE)
head(GTseq_ONLY_popINFO)
dim(GTseq_ONLY_popINFO)

####Input a list of the individuals that have passed all the filters so far:
kept_samples <- read.delim("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/filteredSample_names.txt",header=FALSE)
head(kept_samples)

##############################  RAD LOCI ###################################
############# Make TXT files to split up genepop files by different groups of individuals ####

head(chin_all_popINFO)
dim(chin_all_popINFO)

trimmed_chin_popINFO <- chin_all_popINFO[chin_all_popINFO$SampleName %in% filteredSample_names,]
dim(trimmed_chin_popINFO)

#group A is just the Nushigak
Nushigak <- trimmed_chin_popINFO[trimmed_chin_popINFO$BroadRegion == "Nushigak",1 ]
head(Nushigak)
length(Nushigak)
drop.levels(Nushigak)

#group B is just the Kuskokwim 
Kusksokwim <- trimmed_chin_popINFO[trimmed_chin_popINFO$BroadRegion == "Kusksokwim",1 ]
head(Kusksokwim)
length(Kusksokwim)
drop.levels(Kusksokwim)

#group C is just the Nushigak and the Kuskokwim 
Kuskso_Nush <- trimmed_chin_popINFO[(trimmed_chin_popINFO$BroadRegion == "Kusksokwim" | trimmed_chin_popINFO$BroadRegion == "Nushigak"),1 ]
head(Kuskso_Nush)
length(Kuskso_Nush)
drop.levels(Kuskso_Nush)

#### Write these out to txt files to be used to cut down genepop files 

#A
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/Nushigak.txt", "wb")
write.table(Nushigak,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#B
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/Kusksokwim.txt", "wb")
write.table(Kusksokwim,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#C
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/Kuskso_Nush.txt", "wb")
write.table(Kuskso_Nush,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



