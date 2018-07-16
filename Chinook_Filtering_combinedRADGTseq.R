###  Chinook panel: filtering CombinedRADGTseq data 
###    to keep the loci that are most powerful in a refined panel
###    goal is to differentiate Nushigak and Kusksokwim fish
### Carolyn Tarpey | July 2018 
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

#load the haplotype genepop file that is all the 847 Combined RADGTseq data
chin_comb_RADGTseq <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_comb_RADGTseq)
chin_comb_RADGTseq[1:6,1:6]

#prep for later- the names of the loci and the number: 
chin_comb_RADGTseq_loci <- colnames(chin_comb_RADGTseq)[2:848]
head(chin_comb_RADGTseq_loci)
nloci <- length(chin_comb_RADGTseq_loci)

#a population map that has all the individuals and the populations that they belong to. 
#this file has been edited to combined the Keek pops 
chin_RADGTseq_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/PopulationMap_847RAD_taq_genepop.txt", header =TRUE)
head(chin_RADGTseq_popINFO)
dim(chin_RADGTseq_popINFO)


###########################Combined RAD GTseq starting with 847 loci #####################
####### Merge the Pop info with the genotypes ######
#assign each sample to a population 
chin_comb_RADGTseq_pop<-merge(chin_RADGTseq_popINFO,chin_comb_RADGTseq, by="SampleName")
chin_comb_RADGTseq_pop[1:5,1:15]
dim(chin_comb_RADGTseq_pop)


##############Minor Allele Frequency MAF#############

#Function to get minor allele frequency for each locus
calculateMAF<-function(genotypes){
  genotypeList<-sort(unique(genotypes))
  genotypeList<-genotypeList[genotypeList != "0000"]
  allelesList1<-substr(genotypeList,1,2)
  allelesList2<-substr(genotypeList,3,4)
  allelesList<-unique(c(allelesList1,allelesList2))
  allele1Counts<-sum(str_count(genotypes,allelesList[1]))
  allele2Counts<-sum(str_count(genotypes,allelesList[2]))
  if(length(allelesList)==1){
    MAF=0
  }else if(allele1Counts>=allele2Counts){
    MAF<-allele2Counts/(allele1Counts+allele2Counts)
  }else{
    MAF<-allele1Counts/(allele1Counts+allele2Counts)
  }
  return(MAF)
}



######### Test MAF >= 0.05 in at least one of 17 pops ############
#create a dataframe for our MAF results
npops<- length(unique(chin_comb_RADGTseq_pop$Population))

RADGTseq_POP_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(RADGTseq_POP_MAF)<-as.vector(unique(chin_comb_RADGTseq_pop$Population)) #name the rows by the population names
colnames(RADGTseq_POP_MAF)<-chin_comb_RADGTseq_loci
RADGTseq_POP_MAF[1:6,1:6]
dim(RADGTseq_POP_MAF)


##### Get the MAF for each locus by population 
for (i in 1:17){ #17 is the unique number of populations we have
  tempPop<-(unique(chin_comb_RADGTseq_pop$Population)[i])
  tempGeno<-subset(chin_comb_RADGTseq_pop, Population == tempPop)
  tempGeno[,1:10]
  # dim(tempGeno)
  tempMAF<-apply(tempGeno[,c(6:dim(chin_comb_RADGTseq_pop)[2])],2,calculateMAF)
  RADGTseq_POP_MAF[i,]<-c(tempMAF)
}

RADGTseq_POP_MAF[1:17,1:5]
dim(RADGTseq_POP_MAF)
#colnames(RADGTseq_POP_MAF[,1])

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 16 populations
lociMAF<-vector()
lociMAF<-apply(RADGTseq_POP_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
lociMAF<-data.frame(keynames=names(lociMAF), value=lociMAF, row.names = NULL)
colnames(lociMAF)<-c("Locus", "PopsMAF")
dim(lociMAF)

#if the loci does not pass the test, delete it
locipassedMAF<-vector()
locipassedMAF<-subset(lociMAF, lociMAF[,2]!=0)
dim(locipassedMAF)
head(locipassedMAF)

#what percent of the loci do we KEEP with this filter? 
Percent_locipassedMAF <- ((dim(locipassedMAF)[1])/(dim(lociMAF)[1]))*100
Percent_locipassedMAF

######### Test MAF >= 0.05 in at least one of the 4 main regions: BroadReportingGroup (BroadRegion, column 5) ############
colnames(chin_comb_RADGTseq_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_comb_RADGTseq_pop$Broad_Region))

RADGTseq_BREG_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(RADGTseq_BREG_MAF)<-as.vector(unique(chin_comb_RADGTseq_pop$Broad_Region)) #name the rows by the population names
colnames(RADGTseq_BREG_MAF)<-chin_comb_RADGTseq_loci
dim(RADGTseq_BREG_MAF)
RADGTseq_BREG_MAF[1:3,1:6]

##### Get the MAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_comb_RADGTseq_pop$Broad_Region)[i])
  tempGeno<-subset(chin_comb_RADGTseq_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(6:dim(chin_comb_RADGTseq_pop)[2])],2,calculateMAF)
  RADGTseq_BREG_MAF[i,]<-c(tempMAF)
}

RADGTseq_BREG_MAF[1:3,1:5]
dim(RADGTseq_BREG_MAF)
#colnames(RADGTseq_BREG_MAF[,1])

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_BREGMAF<-vector()
loci_BREGMAF<-apply(RADGTseq_BREG_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
loci_BREGMAF<-data.frame(keynames=names(loci_BREGMAF), value=loci_BREGMAF, row.names = NULL)
colnames(loci_BREGMAF)<-c("Locus", "BroadRegionMAF")
head(loci_BREGMAF)

#if the loci does not pass the test, delete it
locipassed_BREGMAF<-vector()
locipassed_BREGMAF<-subset(loci_BREGMAF, loci_BREGMAF[,2]!=0)
dim(locipassed_BREGMAF)
head(locipassed_BREGMAF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMAF <- ((dim(locipassed_BREGMAF)[1])/847)*100
Percent_locipassed_BREGMAF

#get list of loci that made it past this filter
locipassed_BREGMAF_keepList <- locipassed_BREGMAF$Locus
length(locipassed_BREGMAF_keepList)
head(locipassed_BREGMAF_keepList)
drop.levels(locipassed_BREGMAF_keepList)


################# Filter out Loci that do not meet the MAF threshold of 0.05 in 1 or more regions ##########################
filtered_MAF_Genos<-chin_comb_RADGTseq_pop[,((colnames(chin_comb_RADGTseq_pop) %in% locipassed_BREGMAF_keepList) | (colnames(chin_comb_RADGTseq_pop) %in% colnames(chin_RADGTseq_popINFO)))]
dim(filtered_MAF_Genos)
filtered_MAF_Genos[1:4,1:15]


################## Output lists of Loci that passed MAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/locipassed_BroadREGMAF_keepList.txt", "wb")
write.table(locipassed_BREGMAF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)












