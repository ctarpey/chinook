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

####import a list of the individuals that have passed all the filters so far:
filteredSample_names <- read.table("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/filteredSample_names.txt", header=FALSE)
head(filteredSample_names)
dim(filteredSample_names)

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
RADtaqAmp_Rcombined[846:849,846:849]
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
kept_samples <- read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/filteredSample_names.txt",header=FALSE)
head(kept_samples)

############################## new RAD LOCI ###################################
############# Make TXT files to split up genepop files by different groups of individuals ####
head(chin_RAD_all_popINFO)
dim(chin_RAD_all_popINFO)
length(unique(chin_RAD_all_popINFO$FineReportingGroup))
trimmed_chin_popINFO <- chin_RAD_all_popINFO[chin_RAD_all_popINFO$SampleName %in% filteredSample_names$V1,]
dim(trimmed_chin_popINFO)

#group A is just the Nushigak
Nushigak <- trimmed_chin_popINFO[trimmed_chin_popINFO$BroadRegion == "Nushigak",1 ]
Nushigak <- drop.levels(Nushigak)
head(Nushigak)
length(Nushigak)

#group B is just the Kuskokwim 
Kusksokwim <- trimmed_chin_popINFO[trimmed_chin_popINFO$BroadRegion == "Kusksokwim",1 ]
Kusksokwim <- drop.levels(Kusksokwim)
head(Kusksokwim)
length(Kusksokwim)

#group C is just the Nushigak and the Kuskokwim 
Kuskso_Nush <- trimmed_chin_popINFO[(trimmed_chin_popINFO$BroadRegion == "Kusksokwim" | trimmed_chin_popINFO$BroadRegion == "Nushigak"),1 ]
Kusko_Nush <- drop.levels(Kuskso_Nush)
head(Kuskso_Nush)
length(Kuskso_Nush)

#group D is the all, the with all the populations including yukon, norton, upper kusk, just no cook inlet
NewRAD_AllnoCookpops <- trimmed_chin_popINFO$SampleName
head(NewRAD_AllnoCookpops)

#group E is the populations that are also in the all of the combined RADGTseq, for comparison purposes. 
pops_inRADGT <- unique(RADtaqAmp_popINFO$Population)
pops_inRADGT
head(trimmed_chin_popINFO)
dim(trimmed_chin_popINFO)
newRAD_RADGT_overlap <- trimmed_chin_popINFO[trimmed_chin_popINFO$Population %in% RADtaqAmp_popINFO$Population,] 
dim(newRAD_RADGT_overlap)
newRAD_RADGT_overlap <- newRAD_RADGT_overlap$SampleName

#### Write these out to txt files to be used to cut down genepop files 
#A
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/Nushigak_newRAD.txt", "wb")
write.table(Nushigak,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol=",\n")
close(outputFile)

#B
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/Kusksokwim_newRAD.txt", "wb")
write.table(Kusksokwim,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol=",\n")
close(outputFile)

#C
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/Kusk_Nush_newRAD.txt", "wb")
write.table(Kuskso_Nush,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol=",\n")
close(outputFile)

#D
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/NewRAD_AllnoCookpops.txt", "wb")
write.table(NewRAD_AllnoCookpops,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol=",\n")
close(outputFile)

#E
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/FilteringRawGenotypes/newRAD_RADGT_overlap.txt", "wb")
write.table(newRAD_RADGT_overlap,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol=",\n")
close(outputFile)

##############################  RADGTseq combined LOCI ###################################
############# Make TXT files to split up genepop files by different groups of individuals ####
head(RADtaqAmp_popINFO)
dim(RADtaqAmp_popINFO)
length(unique(RADtaqAmp_popINFO$FineReportingGroup))

#group A is just the Nushigak
Nushigak_radGT <- RADtaqAmp_popINFO[RADtaqAmp_popINFO$Broad_Region == "Nushigak",1 ]
Nushigak_radGT <- drop.levels(Nushigak_radGT)
head(Nushigak_radGT)
length(Nushigak_radGT)

#group B is just the Kuskokwim 
Kusksokwim_radGT <- RADtaqAmp_popINFO[RADtaqAmp_popINFO$Broad_Region == "Kuskokwim",1 ]
Kusksokwim_radGT <- drop.levels(Kusksokwim_radGT)
head(Kusksokwim_radGT)
length(Kusksokwim_radGT)

#group C is just the Nushigak and the Kuskokwim 
Kuskso_Nush_radGT  <- RADtaqAmp_popINFO[(RADtaqAmp_popINFO$Broad_Region == "Kuskokwim" | RADtaqAmp_popINFO$Broad_Region == "Nushigak"),1 ]
Kuskso_Nush_radGT  <- drop.levels(Kuskso_Nush_radGT )
head(Kuskso_Nush_radGT )
length(Kuskso_Nush_radGT )


#### Write these out to txt files to be used to cut down genepop files 

#A
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/Nushigak_radGT.txt", "wb")
write.table(Nushigak_radGT,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#B
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/Kuskokwim_radGT.txt", "wb")
write.table(Kusksokwim_radGT,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#C
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/Kusk_Nush_radGT.txt", "wb")
write.table(Kuskso_Nush_radGT,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)




##############################  GTseqONLY combined LOCI ###################################
############# Make TXT files to split up genepop files by different groups of individuals ####
head(GTseq_ONLY_popINFO)
dim(GTseq_ONLY_popINFO)
length(unique(GTseq_ONLY_popINFO$FineReportingGroup))

#group A is just the Nushigak
Nushigak_GTonly <- GTseq_ONLY_popINFO[GTseq_ONLY_popINFO$Broad_Region == "Nushigak",1 ]
Nushigak_GTonly <- drop.levels(Nushigak_GTonly)
head(Nushigak_GTonly)
length(Nushigak_GTonly)

#group B is just the Kuskokwim 
Kusksokwim_GTonly <- GTseq_ONLY_popINFO[GTseq_ONLY_popINFO$Broad_Region == "Kuskokwim",1 ]
Kusksokwim_GTonly <- drop.levels(Kusksokwim_GTonly)
head(Kusksokwim_GTonly)
length(Kusksokwim_GTonly)

#group C is just the Nushigak and the Kuskokwim 
Kuskso_Nush_GTonly  <- GTseq_ONLY_popINFO[(GTseq_ONLY_popINFO$Broad_Region == "Kuskokwim" | GTseq_ONLY_popINFO$Broad_Region == "Nushigak"),1 ]
Kuskso_Nush_GTonly  <- drop.levels(Kuskso_Nush_GTonly )
head(Kuskso_Nush_GTonly )
length(Kuskso_Nush_GTonly )


#### Write these out to txt files to be used to cut down genepop files 

#A
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/Nushigak_GTonly", "wb")
write.table(Nushigak_GTonly,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#B
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/Kusksokwim_GTonly.txt", "wb")
write.table(Kusksokwim_GTonly,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#C
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/Kuskso_Nush_GTonly.txt", "wb")
write.table(Kuskso_Nush_GTonly,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



