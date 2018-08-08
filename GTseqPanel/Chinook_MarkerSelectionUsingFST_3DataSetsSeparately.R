### To select markers for a refined Chinook panel 
###   Using Ranked FST of one snp per tag loci in a hierarchical structure <--- Method 1
###    
### to identify markers that distinguish the Kusk from Nush
###
###   Carolyn Tarpey | July 2018
### ---------------------------------------

#This code requires the estimated FST generated from Genepop for each of the population groups within each of the data types. 
#The output from genepop was assembled in an excel file and then moved to a text file, a single table with all the FST of each locus, per population group for each of the data types. 

#There are four population groupings per data type:  
#1. All populations
#2. All Kuskokwim populations 
#3. All Nushagak populations
#4. All Kuskokwim populations grouped as one, all Nushigak grouped as the other

## There are three data types
# 1. allCombined are the new RAD loci that have not been in a panel before <- right now using One SNP per tag- the one with the greatest MAF
#    we arent interested in loci that have low MAF, they arent useful on a panel because they are so rare in the natural population. 
# 2. RADGTseq loci are the panel loci that had additional RAD data for other individuals, so they have 1544 individuals
# 3. GTseqONLY are the loci that were on the panel that didnt have additional RAD loci, so they have 995 individuals

#Haplotype and duplicated loci
# the GTseqONLY data have duplicated loci that we have removed from the dataset that we dont know how to get an FST estimate for. 
#the GTseqONLY and the RADGTseq data both have haplotypes mixed in with the other one snp per tag loci. 
#They already have primers made, so they dont have to have the FST per locus 

###### Overlapping loci <----------------------------THIS WAS NOT DONE IN THIS SCRIPT!!!!
##There are no overlapping loci in the GTseqONLY and the combinedRADGTseq, 
##but there is probably some between the combinedRADGTseq and the NewRAD, the allCombined, so they need to be identified
## we will keep the version that is in the combined RADGTseq because they have primers made and incorporate any haplotypes.


### Population Maps 
### This code also requires a population map, which each individual of the analysis in the rows, and their population named in columns. 

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
allCombined_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADallCombined_nonduplicated_FST_results.txt", header=TRUE,  row.names=NULL,
                     stringsAsFactors = FALSE, na.strings = "-" )
allCombined_FST_file[1:5,1:4]
dim(allCombined_FST_file)

## RADGTseq data set estimated FST
RADGTseq_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADGTSeq_nonduplicated_FST_results.txt", header=TRUE,  row.names=NULL,
                     stringsAsFactors = FALSE, na.strings = "-" )
RADGTseq_FST_file[1:5,1:4]
dim(RADGTseq_FST_file)

## GTseqONLY data set estimated FST
GTseqONLY_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/GTSeqONLY_nonduplicated_FST_results.txt", header=TRUE, row.names=NULL,
                     stringsAsFactors = FALSE, na.strings = "-" )
GTseqONLY_FST_file[1:5,1:4]
dim(GTseqONLY_FST_file)


###################################### Marker selection based on ranked FSt of loci 
##look at the FST of the snps and choose the ones that have the highest FST for the primer design and panel development. 

##################################################
############## GTseqONLY_FST_file ################
##################################################

names(GTseqONLY_FST_file)
dim(GTseqONLY_FST_file)
head(GTseqONLY_FST_file)

##########################
##### Create data frames for each of the population groups within GTseqONLY_FST_file data set

##Find the top FST for All
GTseqONLY_FST_ALL <- GTseqONLY_FST_file[,c("row.names","All")]
head(GTseqONLY_FST_ALL)
GTseqONLY_FST_ALL <- GTseqONLY_FST_ALL[order(GTseqONLY_FST_ALL$All, decreasing=TRUE),]
head(GTseqONLY_FST_ALL)

##Find the top FST for the Kusk
GTseqONLY_FST_Kusk <- GTseqONLY_FST_file[,c("row.names","Kusk")]
head(GTseqONLY_FST_Kusk)
GTseqONLY_FST_Kusk <- GTseqONLY_FST_Kusk[order(GTseqONLY_FST_Kusk$Kusk, decreasing=TRUE),]
head(GTseqONLY_FST_Kusk)

##Find the top FST for the Nush
GTseqONLY_FST_Nush <- GTseqONLY_FST_file[,c("row.names","Nush")]
head(GTseqONLY_FST_Nush)
GTseqONLY_FST_Nush <- GTseqONLY_FST_Nush[order(GTseqONLY_FST_Nush$Nush, decreasing=TRUE),]
head(GTseqONLY_FST_Nush)

##Find the top FST for the Kusk_Nush
GTseqONLY_FST_Kusk_Nush <- GTseqONLY_FST_file[,c("row.names","Kusk_Nush")]
head(GTseqONLY_FST_Kusk_Nush)
GTseqONLY_FST_Kusk_Nush <- GTseqONLY_FST_Kusk_Nush[order(GTseqONLY_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
head(GTseqONLY_FST_Kusk_Nush)


##################################################
############## RADGTseq_FST_file #################
##################################################

names(RADGTseq_FST_file)
dim(RADGTseq_FST_file)
head(RADGTseq_FST_file)

##########################
##### Create data frames for each of the population groups within RADGTseq_FST_file data set

##Find the top FST for All
RADGTseq_FST_ALL <- RADGTseq_FST_file[,c("row.names","All")]
head(RADGTseq_FST_ALL)
RADGTseq_FST_ALL <- RADGTseq_FST_ALL[order(RADGTseq_FST_ALL$All, decreasing=TRUE),]
head(RADGTseq_FST_ALL)

##Find the top FST for the Kusk
RADGTseq_FST_Kusk <- RADGTseq_FST_file[,c("row.names","Kusk")]
head(RADGTseq_FST_Kusk)
RADGTseq_FST_Kusk <- RADGTseq_FST_Kusk[order(RADGTseq_FST_Kusk$Kusk, decreasing=TRUE),]
head(RADGTseq_FST_Kusk)

##Find the top FST for the Nush
RADGTseq_FST_Nush <- RADGTseq_FST_file[,c("row.names","Nush")]
head(RADGTseq_FST_Nush)
RADGTseq_FST_Nush <- RADGTseq_FST_Nush[order(RADGTseq_FST_Nush$Nush, decreasing=TRUE),]
head(RADGTseq_FST_Nush)

##Find the top FST for the Kusk_Nush
RADGTseq_FST_Kusk_Nush <- RADGTseq_FST_file[,c("row.names","Kusk_Nush")]
head(RADGTseq_FST_Kusk_Nush)
RADGTseq_FST_Kusk_Nush <- RADGTseq_FST_Kusk_Nush[order(RADGTseq_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
head(RADGTseq_FST_Kusk_Nush)

##################################################
############## allCombined_FST_file ############## 
##################################################

names(allCombined_FST_file)
dim(allCombined_FST_file)
head(allCombined_FST_file)

##########################
##### Create data frames for each of the population groups within allCombined_FST_file data set

##Find the top FST for All
allCombined_FST_ALL <- allCombined_FST_file[,c("row.names","All")]
head(allCombined_FST_ALL)
allCombined_FST_ALL <- allCombined_FST_ALL[order(allCombined_FST_ALL$All, decreasing=TRUE),]
head(allCombined_FST_ALL)

##Find the top FST for the Kusk
allCombined_FST_Kusk <- allCombined_FST_file[,c("row.names","Kusk")]
head(allCombined_FST_Kusk)
allCombined_FST_Kusk <- allCombined_FST_Kusk[order(allCombined_FST_Kusk$Kusk, decreasing=TRUE),]
head(allCombined_FST_Kusk)

##Find the top FST for the Nush
allCombined_FST_Nush <- allCombined_FST_file[,c("row.names","Nush")]
head(allCombined_FST_Nush)
allCombined_FST_Nush <- allCombined_FST_Nush[order(allCombined_FST_Nush$Nush, decreasing=TRUE),]
head(allCombined_FST_Nush)

##Find the top FST for the Kusk_Nush
allCombined_FST_Kusk_Nush <- allCombined_FST_file[,c("row.names","Kusk_Nush")]
head(allCombined_FST_Kusk_Nush)
allCombined_FST_Kusk_Nush <- allCombined_FST_Kusk_Nush[order(allCombined_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
head(allCombined_FST_Kusk_Nush)


########################################################################
##plot the range of the FST for each of the groups

##################################################
############## GTseqONLY_FST_file ################
##################################################

GTseqONLY_FST_ALL_ranked <- GTseqONLY_FST_ALL[order(GTseqONLY_FST_ALL$All, decreasing=TRUE),]
GTseqONLY_FST_ALL_ranked$rank<-seq(1,dim(GTseqONLY_FST_ALL_ranked)[1],by=1)
ggplot()+geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All),color ="goldenrod1" )+ggtitle("FST range of GTSeqONLY ALL populations ")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

GTseqONLY_FST_Kusk_ranked <- GTseqONLY_FST_Kusk[order(GTseqONLY_FST_Kusk$Kusk, decreasing=TRUE),]
GTseqONLY_FST_Kusk_ranked$rank<-seq(1,dim(GTseqONLY_FST_Kusk_ranked)[1],by=1)
ggplot()+geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk),color ="violetred1" )+ggtitle("FST range of GTSeqONLY Kuskokwim populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

GTseqONLY_FST_Nush_ranked <- GTseqONLY_FST_Nush[order(GTseqONLY_FST_Nush$Nush, decreasing=TRUE),]
GTseqONLY_FST_Nush_ranked$rank<-seq(1,dim(GTseqONLY_FST_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush),color ="orangered1" )+ggtitle("FST range of GTSeqONLY Nushagak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

GTseqONLY_FST_Kusk_Nush_ranked <- GTseqONLY_FST_Kusk_Nush[order(GTseqONLY_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
GTseqONLY_FST_Kusk_Nush_ranked$rank<-seq(1,dim(GTseqONLY_FST_Kusk_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush),color ="red1" )+ggtitle("FST range of GTSeqONLY Kuskokwim and Nushigak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

#plot all four on the same graph in different colors
ggplot()+geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All, colour = "GTseqONLY_ALL")) +
  geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "GTseqONLY_Kusk")) +
  geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "GTseqONLY_Nush")) +
  geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "GTseqONLY_Kusk_Nush")) +
  ggtitle("FST range of the four GTSeqONLY population groups") +theme_bw() + 
  scale_colour_manual(values = c("goldenrod1", "violetred1", "orangered1", "red1"), 
                      breaks = c("GTseqONLY_ALL", "GTseqONLY_Kusk", "GTseqONLY_Nush", "GTseqONLY_Kusk_Nush")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for GTseqONLY population groups", x = "Rank", y = "FST", size = 20))  +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("goldenrod1", "violetred1", "orangered1", "red1"))))

##################################################
############## RADGTseq_FST_file #################
##################################################

RADGTseq_FST_ALL_ranked <- RADGTseq_FST_ALL[order(RADGTseq_FST_ALL$All, decreasing=TRUE),]
RADGTseq_FST_ALL_ranked$rank<-seq(1,dim(RADGTseq_FST_ALL_ranked)[1],by=1)
ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All),color ="slateblue1" )+ggtitle("FST range of RADGTseq ALL populations ")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

RADGTseq_FST_Kusk_ranked <- RADGTseq_FST_Kusk[order(RADGTseq_FST_Kusk$Kusk, decreasing=TRUE),]
RADGTseq_FST_Kusk_ranked$rank<-seq(1,dim(RADGTseq_FST_Kusk_ranked)[1],by=1)
ggplot()+geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk),color ="blue")+ggtitle("FST range of RADGTseq Kuskokwim populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

RADGTseq_FST_Nush_ranked <- RADGTseq_FST_Nush[order(RADGTseq_FST_Nush$Nush, decreasing=TRUE),]
RADGTseq_FST_Nush_ranked$rank<-seq(1,dim(RADGTseq_FST_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush),color ="deepskyblue2")+ggtitle("FST range of RADGTseq Nushagak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

RADGTseq_FST_Kusk_Nush_ranked <- RADGTseq_FST_Kusk_Nush[order(RADGTseq_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
RADGTseq_FST_Kusk_Nush_ranked$rank<-seq(1,dim(RADGTseq_FST_Kusk_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush),color ="royalblue2")+ggtitle("FST range of RADGTseq Kuskokwim and Nushigak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")


#plot all four on the same graph in different colors
ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All, colour = "RADGTseq_ALL")) +
  geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "RADGTseq_Kusk")) +
  geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "RADGTseq_Nush")) +
  geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "RADGTseq_Kusk_Nush")) + theme_bw() +
  scale_colour_manual(values = c("slateblue1", "blue", "deepskyblue2", "royalblue2"), 
                      breaks = c("RADGTseq_ALL","RADGTseq_Kusk","RADGTseq_Nush","RADGTseq_Kusk_Nush")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for RADGTseq population groups", x = "Rank", y = "FST", size = 20)) + 
  
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("slateblue1", "blue", "deepskyblue2", "royalblue2"))))



##################################################
############## allCombined_FST_file ############## 
##################################################

allCombined_FST_ALL_ranked <- allCombined_FST_ALL[order(allCombined_FST_ALL$All, decreasing=TRUE),]
allCombined_FST_ALL_ranked$rank<-seq(1,dim(allCombined_FST_ALL_ranked)[1],by=1)
ggplot()+geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All), color="darkgreen" )+ggtitle("FST range of allCombined ALL populations ")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

allCombined_FST_Kusk_ranked <- allCombined_FST_Kusk[order(allCombined_FST_Kusk$Kusk, decreasing=TRUE),]
allCombined_FST_Kusk_ranked$rank<-seq(1,dim(allCombined_FST_Kusk_ranked)[1],by=1)
ggplot()+geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk), color="green4" )+ggtitle("FST range of allCombined Kuskokwim populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

allCombined_FST_Nush_ranked <- allCombined_FST_Nush[order(allCombined_FST_Nush$Nush, decreasing=TRUE),]
allCombined_FST_Nush_ranked$rank<-seq(1,dim(allCombined_FST_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush), color="olivedrab4" )+ggtitle("FST range of allCombined Nushagak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

allCombined_FST_Kusk_Nush_ranked <- allCombined_FST_Kusk_Nush[order(allCombined_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
allCombined_FST_Kusk_Nush_ranked$rank<-seq(1,dim(allCombined_FST_Kusk_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush), color="springgreen" )+ggtitle("FST range of allCombined Kuskokwim and Nushigak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

#plot all four on the same graph in different colors
ggplot()+geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All, colour= "allCombined_ALL")) +
  geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour= "allCombined_Kusk")) +
  geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush, colour= "allCombined_Nush")) +
  geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour= "allCombined_Kusk_Nush")) + theme_bw() + 
  scale_colour_manual(values = c("darkgreen", "green4", "olivedrab4", "springgreen"), 
                      breaks = c("allCombined_ALL","allCombined_Kusk","allCombined_Nush","allCombined_Kusk_Nush")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for allCombined population groups", x = "Rank", y = "FST", size = 20)) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("darkgreen", "green4", "olivedrab4", "springgreen"))))

########################################################################

##plot all of the FST ranks on the same plot

ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All, colour = "RADGTseq_ALL"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "RADGTseq_Kusk"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "RADGTseq_Nush"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "RADGTseq_Kusk_Nush"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All, colour = "allCombined_ALL"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "allCombined_Kusk"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "allCombined_Nush"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "allCombined_Kusk_Nush"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All, colour = "GTseqONLY_ALL"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "GTseqONLY_Kusk"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "GTseqONLY_Nush"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "GTseqONLY_Kusk_Nush"), size = 1, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() +
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "darkgreen", "green4", "olivedrab4", "springgreen", "goldenrod1", "violetred1", "orangered1", "red1","slateblue1", "blue", "deepskyblue2", "royalblue2") ,
  breaks = c("RADGTseq_ALL", "RADGTseq_Kusk","RADGTseq_Nush","RADGTseq_Kusk_Nush","GTseqONLY_ALL","GTseqONLY_Kusk", "GTseqONLY_Nush","GTseqONLY_Kusk_Nush",
             "allCombined_ALL","allCombined_Kusk", "allCombined_Nush","allCombined_Kusk_Nush" )) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("slateblue1", "blue", "deepskyblue2", "royalblue2", "goldenrod1", "violetred1", "orangered1", "red1","darkgreen", "green4", "olivedrab4", "springgreen"))))



##zoom in + xlim(0, 1000) + ylim(-0.01, 0.80) + 
ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All, colour = "RADGTseq_ALL"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "RADGTseq_Kusk"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "RADGTseq_Nush"), size = 1, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "RADGTseq_Kusk_Nush"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All, colour = "allCombined_ALL"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "allCombined_Kusk"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "allCombined_Nush"), size = 1, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "allCombined_Kusk_Nush"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All, colour = "GTseqONLY_ALL"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "GTseqONLY_Kusk"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "GTseqONLY_Nush"), size = 1, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "GTseqONLY_Kusk_Nush"), size = 1, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() + xlim(0, 1000) + ylim(-0.01, 0.80) + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "darkgreen", "green4", "olivedrab4", "springgreen", "goldenrod1", "violetred1", "orangered1", "red1","slateblue1", "blue", "deepskyblue2", "royalblue2") ,
                      breaks = c("RADGTseq_ALL", "RADGTseq_Kusk","RADGTseq_Nush","RADGTseq_Kusk_Nush","GTseqONLY_ALL","GTseqONLY_Kusk", "GTseqONLY_Nush","GTseqONLY_Kusk_Nush",
                                 "allCombined_ALL","allCombined_Kusk", "allCombined_Nush","allCombined_Kusk_Nush" )) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("slateblue1", "blue", "deepskyblue2", "royalblue2", "goldenrod1", "violetred1", "orangered1", "red1","darkgreen", "green4", "olivedrab4", "springgreen"))))


##zoom in even more 
ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All, colour = "RADGTseq_ALL"), size = 2, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "RADGTseq_Kusk"), size = 2, shape = 20) +
  geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "RADGTseq_Nush"), size = 2, shape = 20) +
  geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "RADGTseq_Kusk_Nush"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All, colour = "allCombined_ALL"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "allCombined_Kusk"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "allCombined_Nush"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "allCombined_Kusk_Nush"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All, colour = "GTseqONLY_ALL"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "GTseqONLY_Kusk"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "GTseqONLY_Nush"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "GTseqONLY_Kusk_Nush"), size = 2, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() + xlim(0, 250) + ylim(-0.01, 0.50) + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "darkgreen", "green4", "olivedrab4", "springgreen", "goldenrod1", "violetred1", "orangered1", "red1","slateblue1", "blue", "deepskyblue2", "royalblue2") ,
                      breaks = c("RADGTseq_ALL", "RADGTseq_Kusk","RADGTseq_Nush","RADGTseq_Kusk_Nush","GTseqONLY_ALL","GTseqONLY_Kusk", "GTseqONLY_Nush","GTseqONLY_Kusk_Nush",
                                 "allCombined_ALL","allCombined_Kusk", "allCombined_Nush","allCombined_Kusk_Nush" )) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("slateblue1", "blue", "deepskyblue2", "royalblue2", "goldenrod1", "violetred1", "orangered1", "red1","darkgreen", "green4", "olivedrab4", "springgreen"))))


########################################################################

##plot all of the FST ranks on the same plot####<------------------------------------------------START HERE, do the sets individually

ggplot()+geom_point(data=RADGTseq_FST_ALL_ranked,aes(x=rank,y=All, colour = "RADGTseq_ALL"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_ALL_ranked,aes(x=rank,y=All, colour = "allCombined_ALL"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_ALL_ranked,aes(x=rank,y=All, colour = "GTseqONLY_ALL"), size = 2, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() +
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "darkgreen",  "goldenrod1", "slateblue1" ) ,
    breaks = c("RADGTseq_ALL", "GTseqONLY_ALL", "allCombined_ALL")) + xlim(0, 250) + ylim(-0.01, 0.50) + 
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("slateblue1","goldenrod1", "darkgreen"))))



##zoom in + xlim(0, 1000) + ylim(-0.01, 0.80) + 
ggplot() +
  geom_point(data=RADGTseq_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "RADGTseq_Kusk"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "allCombined_Kusk"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_ranked,aes(x=rank,y=Kusk, colour = "GTseqONLY_Kusk"), size = 2, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() + xlim(0, 250) + ylim(-0.01, 0.50) + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all Kuskokwim population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "green4", "violetred1", "blue") ,
                      breaks = c("RADGTseq_Kusk","GTseqONLY_Kusk","allCombined_Kusk")) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("blue", "violetred1","green4"))))


##zoom in even more 
ggplot()+
  geom_point(data=RADGTseq_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "RADGTseq_Nush"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "allCombined_Nush"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Nush_ranked,aes(x=rank,y=Nush, colour = "GTseqONLY_Nush"), size = 2, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() + xlim(0, 250) + ylim(-0.01, 0.50) + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all Nushagak population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "olivedrab4", "orangered1", "deepskyblue2") ,
                      breaks = c("RADGTseq_Nush", "GTseqONLY_Nush", "allCombined_Nush")) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("deepskyblue2", "orangered1","olivedrab4"))))




##zoom in even more 
ggplot()+
  geom_point(data=RADGTseq_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "RADGTseq_Kusk_Nush"), size = 2, shape = 20) +
  geom_point(data=allCombined_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "allCombined_Kusk_Nush"), size = 2, shape = 20) +
  geom_point(data=GTseqONLY_FST_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "GTseqONLY_Kusk_Nush"), size = 2, shape = 20) +
  geom_hline(yintercept = 0, color = "black", size= .5) + theme_bw() + xlim(0, 250) + ylim(-0.01, 0.50) + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range for all Nushagak and Kuskokwim population groups", x = "Rank", y = "FST", size = 20)) + 
  scale_colour_manual(values =  c( "springgreen","red1", "royalblue2") ,
                      breaks = c("RADGTseq_Kusk_Nush","GTseqONLY_Kusk_Nush","allCombined_Kusk_Nush" )) +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c( "royalblue2", "red1", "springgreen"))))













##################################################
############## GTseqONLY_FST_file ################
##################################################

################# Choose the top FST loci in each of the population groups for the data set 
#top FST loci for GTseqONLY_FST_ALL
#change the variable to the number you want 
GTseqONLY_FST_ALL_var <- 25 #<--------------------------------------------dictates the number of loci
GTseqONLY_FST_ALL_set <- GTseqONLY_FST_ALL[1:GTseqONLY_FST_ALL_var,]
GTseqONLY_FST_ALL_set_avgFST <- mean(GTseqONLY_FST_ALL_set$All)
GTseqONLY_FST_ALL_set_maxFST <- max(GTseqONLY_FST_ALL_set$All)
GTseqONLY_FST_ALL_set_minFST <- min(GTseqONLY_FST_ALL_set$All)
cat("GTseqONLY_FST_ALL FST max, min, average: ",GTseqONLY_FST_ALL_set_maxFST, GTseqONLY_FST_ALL_set_minFST, GTseqONLY_FST_ALL_set_avgFST )

#top FST loci for GTseqONLY_FST_Kusk
GTseqONLY_FST_Kusk_var <- 25 #<--------------------------------------------dictates the number of loci
GTseqONLY_FST_Kusk_set <- GTseqONLY_FST_Kusk[1:GTseqONLY_FST_Kusk_var,]
GTseqONLY_FST_Kusk_set_avgFST <- mean(GTseqONLY_FST_Kusk_set$Kusk)
GTseqONLY_FST_Kusk_set_maxFST <- max(GTseqONLY_FST_Kusk_set$Kusk)
GTseqONLY_FST_Kusk_set_minFST <- min(GTseqONLY_FST_Kusk_set$Kusk)
cat("GTseqONLY_FST_Kusk FST max, min, average: ", GTseqONLY_FST_Kusk_set_maxFST, GTseqONLY_FST_Kusk_set_minFST, GTseqONLY_FST_Kusk_set_avgFST )

#top FST loci for GTseqONLY_FST_Nush
#change the variable to the number you want 
GTseqONLY_FST_Nush_var <- 25 #<--------------------------------------------dictates the number of loci
GTseqONLY_FST_Nush_set <- GTseqONLY_FST_Nush[1:GTseqONLY_FST_Nush_var,]
GTseqONLY_FST_Nush_set_avgFST <- mean(GTseqONLY_FST_Nush_set$Nush)
GTseqONLY_FST_Nush_set_maxFST <- max(GTseqONLY_FST_Nush_set$Nush)
GTseqONLY_FST_Nush_set_minFST <- min(GTseqONLY_FST_Nush_set$Nush)
cat("GTseqONLY_FST_Nush FST max, min, average: ",GTseqONLY_FST_Nush_set_maxFST, GTseqONLY_FST_Nush_set_minFST, GTseqONLY_FST_Nush_set_avgFST )

#top FST loci for GTseqONLY_FST_Kusk_Nush
GTseqONLY_FST_Kusk_Nush_var <- 25 #<--------------------------------------------dictates the number of loci
GTseqONLY_FST_Kusk_Nush_set <- GTseqONLY_FST_Kusk_Nush[1:GTseqONLY_FST_Kusk_Nush_var,]
GTseqONLY_FST_Kusk_Nush_set_avgFST <- mean(GTseqONLY_FST_Kusk_Nush_set$Kusk_Nush)
GTseqONLY_FST_Kusk_Nush_set_maxFST <- max(GTseqONLY_FST_Kusk_Nush_set$Kusk_Nush)
GTseqONLY_FST_Kusk_Nush_set_minFST <- min(GTseqONLY_FST_Kusk_Nush_set$Kusk_Nush)
cat("GTseqONLY_FST_Kusk_Nush FST max, min, average: ", GTseqONLY_FST_Kusk_Nush_set_maxFST, GTseqONLY_FST_Kusk_Nush_set_minFST, GTseqONLY_FST_Kusk_Nush_set_avgFST )

### see which loci are in the four of the population groups
GTseqONLY_FST_union <- c(GTseqONLY_FST_Kusk_Nush_set$row.names, GTseqONLY_FST_Nush_set$row.names, GTseqONLY_FST_ALL_set$row.names, GTseqONLY_FST_Kusk_set$row.names)
GTseqONLY_FST_union <- unique(GTseqONLY_FST_union)
length(GTseqONLY_FST_union)

###set operations to see which ones overlap, you can only do two sets at a time.  
GTseqONLY_FST_inter_1 <- intersect(GTseqONLY_FST_Kusk_Nush_set$row.names, GTseqONLY_FST_Nush_set$row.names)
GTseqONLY_FST_inter_2 <- intersect(GTseqONLY_FST_inter_1, GTseqONLY_FST_ALL_set$row.names)
GTseqONLY_FST_inter <- intersect(GTseqONLY_FST_inter_2, GTseqONLY_FST_Kusk_set$row.names)
length(GTseqONLY_FST_inter_1)
length(GTseqONLY_FST_inter_2)
length(GTseqONLY_FST_inter)


## make a new dataframe that has the FST info for the union of the sets 
GTseqONLY_FST_union_FST <- GTseqONLY_FST_file[GTseqONLY_FST_file$row.names %in% GTseqONLY_FST_union,]  
GTseqONLY_FST_union_FST_sort <-GTseqONLY_FST_union_FST[order(GTseqONLY_FST_union_FST$All, decreasing = TRUE),]
head( GTseqONLY_FST_union_FST_sort)
dim(GTseqONLY_FST_union_FST_sort)

### Write the FST table for the union of the GTseqONLY_FST_union_FST to a file 
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/GTseqONLY_FST_union.txt", "wb")
write.table(GTseqONLY_FST_union_FST,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

##################################################
############## RADGTseq_FST_file #################
##################################################

################# Choose the top FST loci in each of the population groups for the data set 
#top FST loci forRADGTseq_FST_ALL
#change the variable to the number you want 
RADGTseq_FST_ALL_var <- 50 #<--------------------------------------------dictates the number of loci
RADGTseq_FST_ALL_set <- RADGTseq_FST_ALL[1:RADGTseq_FST_ALL_var,]
RADGTseq_FST_ALL_set_avgFST <- mean(RADGTseq_FST_ALL_set$All)
RADGTseq_FST_ALL_set_maxFST <- max(RADGTseq_FST_ALL_set$All)
RADGTseq_FST_ALL_set_minFST <- min(RADGTseq_FST_ALL_set$All)
cat("RADGTseq_FST_ALL FST max, min, average: ", RADGTseq_FST_ALL_set_maxFST, RADGTseq_FST_ALL_set_minFST, RADGTseq_FST_ALL_set_avgFST )

#top FST loci for RADGTseq_FST_Kusk
RADGTseq_FST_Kusk_var <- 50 #<--------------------------------------------dictates the number of loci
RADGTseq_FST_Kusk_set <- RADGTseq_FST_Kusk[1:RADGTseq_FST_Kusk_var,]
RADGTseq_FST_Kusk_set_avgFST <- mean(RADGTseq_FST_Kusk_set$Kusk)
RADGTseq_FST_Kusk_set_maxFST <- max(RADGTseq_FST_Kusk_set$Kusk)
RADGTseq_FST_Kusk_set_minFST <- min(RADGTseq_FST_Kusk_set$Kusk)
cat("RADGTseq_FST_Kusk FST max, min, average: ", RADGTseq_FST_Kusk_set_maxFST, RADGTseq_FST_Kusk_set_minFST, RADGTseq_FST_Kusk_set_avgFST )

#top FST loci for RADGTseq_FST_Nush
#change the variable to the number you want 
RADGTseq_FST_Nush_var <- 50 #<--------------------------------------------dictates the number of loci
RADGTseq_FST_Nush_set <- RADGTseq_FST_Nush[1:RADGTseq_FST_Nush_var,]
RADGTseq_FST_Nush_set_avgFST <- mean(RADGTseq_FST_Nush_set$Nush)
RADGTseq_FST_Nush_set_maxFST <- max(RADGTseq_FST_Nush_set$Nush)
RADGTseq_FST_Nush_set_minFST <- min(RADGTseq_FST_Nush_set$Nush)
cat("RADGTseq_FST_Nush FST max, min, average: ",RADGTseq_FST_Nush_set_maxFST, RADGTseq_FST_Nush_set_minFST, RADGTseq_FST_Nush_set_avgFST )

#top FST loci for RADGTseq_FST_Kusk_Nush
RADGTseq_FST_Kusk_Nush_var <- 150 #<--------------------------------------------dictates the number of loci
RADGTseq_FST_Kusk_Nush_set <- RADGTseq_FST_Kusk_Nush[1:RADGTseq_FST_Kusk_Nush_var,]
RADGTseq_FST_Kusk_Nush_set_avgFST <- mean(RADGTseq_FST_Kusk_Nush_set$Kusk_Nush)
RADGTseq_FST_Kusk_Nush_set_maxFST <- max(RADGTseq_FST_Kusk_Nush_set$Kusk_Nush)
RADGTseq_FST_Kusk_Nush_set_minFST <- min(RADGTseq_FST_Kusk_Nush_set$Kusk_Nush)
cat("RADGTseq_FST_Kusk_Nush FST max, min, average: ", RADGTseq_FST_Kusk_Nush_set_maxFST, RADGTseq_FST_Kusk_Nush_set_minFST, RADGTseq_FST_Kusk_Nush_set_avgFST )

### see which loci are in the four of the population groups
RADGTseq_FST_union <- c(RADGTseq_FST_Kusk_Nush_set$row.names, RADGTseq_FST_Nush_set$row.names, RADGTseq_FST_ALL_set$row.names, RADGTseq_FST_Kusk_set$row.names)
RADGTseq_FST_union <- unique(RADGTseq_FST_union)
length(RADGTseq_FST_union)

###set operations to see which ones overlap, you can only do two sets at a time.  
RADGTseq_FST_inter_1 <- intersect(RADGTseq_FST_Kusk_Nush_set$row.names, RADGTseq_FST_Nush_set$row.names)
RADGTseq_FST_inter_2 <- intersect(RADGTseq_FST_inter_1, RADGTseq_FST_ALL_set$row.names)
RADGTseq_FST_inter <- intersect(RADGTseq_FST_inter_2, RADGTseq_FST_Kusk_set$row.names)
length(RADGTseq_FST_inter_1)
length(RADGTseq_FST_inter_2)
length(RADGTseq_FST_inter)


## make a new dataframe that has the FST info for the union of the sets 
RADGTseq_FST_union_FST <- RADGTseq_FST_file[RADGTseq_FST_file$row.names %in% RADGTseq_FST_union,]  
RADGTseq_FST_union_FST_sort <-RADGTseq_FST_union_FST[order(RADGTseq_FST_union_FST$All, decreasing = TRUE),]
head(RADGTseq_FST_union_FST_sort)
dim(RADGTseq_FST_union_FST_sort)

### Write the FST table for the union of the RADGTseq_FST_union_FST to a file 
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADGTseq_FST_union.txt", "wb")
write.table(RADGTseq_FST_union_FST,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

##################################################
############## allCombined_FST_file ############## 
##################################################

names(allCombined_FST_file)
dim(allCombined_FST_file)
head(allCombined_FST_file)

##########################
##### Create data frames for each of the population groups within allCombined_FST_file data set

##Find the top FST for All
allCombined_FST_ALL <- allCombined_FST_file[,c("row.names","All")]
head(allCombined_FST_ALL)
allCombined_FST_ALL <- allCombined_FST_ALL[order(allCombined_FST_ALL$All, decreasing=TRUE),]
head(allCombined_FST_ALL)

##Find the top FST for the Kusk
allCombined_FST_Kusk <- allCombined_FST_file[,c("row.names","Kusk")]
head(allCombined_FST_Kusk)
allCombined_FST_Kusk <- allCombined_FST_Kusk[order(allCombined_FST_Kusk$Kusk, decreasing=TRUE),]
head(allCombined_FST_Kusk)

##Find the top FST for the Nush
allCombined_FST_Nush <- allCombined_FST_file[,c("row.names","Nush")]
head(allCombined_FST_Nush)
allCombined_FST_Nush <- allCombined_FST_Nush[order(allCombined_FST_Nush$Nush, decreasing=TRUE),]
head(allCombined_FST_Nush)

##Find the top FST for the Kusk_Nush
allCombined_FST_Kusk_Nush <- allCombined_FST_file[,c("row.names","Kusk_Nush")]
head(allCombined_FST_Kusk_Nush)
allCombined_FST_Kusk_Nush <- allCombined_FST_Kusk_Nush[order(allCombined_FST_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
head(allCombined_FST_Kusk_Nush)



################# Choose the top FST loci in each of the population groups for the data set 
#top FST loci forallCombined_FST_ALL
#change the variable to the number you want 
allCombined_FST_ALL_var <- 50 #<--------------------------------------------dictates the number of loci
allCombined_FST_ALL_set <- allCombined_FST_ALL[1:allCombined_FST_ALL_var,]
allCombined_FST_ALL_set_avgFST <- mean(allCombined_FST_ALL_set$All)
allCombined_FST_ALL_set_maxFST <- max(allCombined_FST_ALL_set$All)
allCombined_FST_ALL_set_minFST <- min(allCombined_FST_ALL_set$All)
cat("allCombined_FST_ALL FST max, min, average: ", allCombined_FST_ALL_set_maxFST, allCombined_FST_ALL_set_minFST, allCombined_FST_ALL_set_avgFST )

#top FST loci for allCombined_FST_Kusk
allCombined_FST_Kusk_var <- 50 #<--------------------------------------------dictates the number of loci
allCombined_FST_Kusk_set <- allCombined_FST_Kusk[1:allCombined_FST_Kusk_var,]
allCombined_FST_Kusk_set_avgFST <- mean(allCombined_FST_Kusk_set$Kusk)
allCombined_FST_Kusk_set_maxFST <- max(allCombined_FST_Kusk_set$Kusk)
allCombined_FST_Kusk_set_minFST <- min(allCombined_FST_Kusk_set$Kusk)
cat("allCombined_FST_Kusk FST max, min, average: ", allCombined_FST_Kusk_set_maxFST, allCombined_FST_Kusk_set_minFST, allCombined_FST_Kusk_set_avgFST )

#top FST loci for allCombined_FST_Nush
#change the variable to the number you want 
allCombined_FST_Nush_var <- 50 #<--------------------------------------------dictates the number of loci
allCombined_FST_Nush_set <- allCombined_FST_Nush[1:allCombined_FST_Nush_var,]
allCombined_FST_Nush_set_avgFST <- mean(allCombined_FST_Nush_set$Nush)
allCombined_FST_Nush_set_maxFST <- max(allCombined_FST_Nush_set$Nush)
allCombined_FST_Nush_set_minFST <- min(allCombined_FST_Nush_set$Nush)
cat("allCombined_FST_Nush FST max, min, average: ",allCombined_FST_Nush_set_maxFST, allCombined_FST_Nush_set_minFST, allCombined_FST_Nush_set_avgFST )

#top FST loci for allCombined_FST_Kusk_Nush
allCombined_FST_Kusk_Nush_var <- 150 #<--------------------------------------------dictates the number of loci
allCombined_FST_Kusk_Nush_set <- allCombined_FST_Kusk_Nush[1:allCombined_FST_Kusk_Nush_var,]
allCombined_FST_Kusk_Nush_set_avgFST <- mean(allCombined_FST_Kusk_Nush_set$Kusk_Nush)
allCombined_FST_Kusk_Nush_set_maxFST <- max(allCombined_FST_Kusk_Nush_set$Kusk_Nush)
allCombined_FST_Kusk_Nush_set_minFST <- min(allCombined_FST_Kusk_Nush_set$Kusk_Nush)
cat("allCombined_FST_Kusk_Nush FST max, min, average: ", allCombined_FST_Kusk_Nush_set_maxFST, allCombined_FST_Kusk_Nush_set_minFST, allCombined_FST_Kusk_Nush_set_avgFST )

### see which loci are in the four of the population groups
allCombined_FST_union <- c(allCombined_FST_Kusk_Nush_set$row.names, allCombined_FST_Nush_set$row.names, allCombined_FST_ALL_set$row.names, allCombined_FST_Kusk_set$row.names)
allCombined_FST_union <- unique(allCombined_FST_union)
length(allCombined_FST_union)

###set operations to see which ones overlap, you can only do two sets at a time.  
allCombined_FST_incommon <- c(allCombined_FST_Kusk_Nush_set$row.names, allCombined_FST_Nush_set$row.names, allCombined_FST_ALL_set$row.names, allCombined_FST_Kusk_set$row.names)
allCombined_FST_incommon <- allCombined_FST_incommon[duplicated(allCombined_FST_incommon)]
length(allCombined_FST_incommon)

allCombined_FST_inter_1 <- intersect(allCombined_FST_Kusk_Nush_set$row.names, allCombined_FST_Nush_set$row.names)
allCombined_FST_inter_2 <- intersect(allCombined_FST_inter_1, allCombined_FST_ALL_set$row.names)
allCombined_FST_inter <- intersect(allCombined_FST_inter_2, allCombined_FST_Kusk_set$row.names)
length(allCombined_FST_inter_1)
length(allCombined_FST_inter_2)
length(allCombined_FST_inter)

########################################################
########################################################
########################################################
##I have not done the below ,but it is useful if we want to output some tables of these loci. 

## make a new dataframe that has the FST info for the union of the sets 
allCombined_FST_union_FST <- allCombined_FST_file[allCombined_FST_file$row.names %in% allCombined_FST_union,]  
allCombined_FST_union_FST_sort <-allCombined_FST_union_FST[order(allCombined_FST_union_FST$All, decreasing = TRUE),]
head(allCombined_FST_union_FST_sort)
dim(allCombined_FST_union_FST_sort)

### Write the FST table for the union of the allCombined_FST_union_FST to a file 
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/allCombined_FST_union.txt", "wb")
write.table(allCombined_FST_union_FST,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)



## Change the name of the file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_failed_passed_list.txt", "wb")
write.table(Even_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

Odd_failed_passed_list <- append(NA_ASIA_odd_union_FST_sort$Locus, Odd_failed_union)
head(Odd_failed_passed_list)
length(Odd_failed_passed_list)
## Change the name of the file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_failed_passed_list.txt", "wb")
write.table(Odd_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)
