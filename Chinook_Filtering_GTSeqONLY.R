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

#load the haplotype genepop file that is all the 847 Combined RADGTseq data
chin_GTseq_ONLY <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/combinedPrimerProbe_amplicon_polyResults_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_GTseq_ONLY)
chin_GTseq_ONLY[1:6,1:6]

#Transpose the data so that each individual is a row
columns <- chin_GTseq_ONLY[,1]
SampleName <- colnames(chin_GTseq_ONLY)[-1]
chin_GTseq_ONLY <- as.data.frame(t(chin_GTseq_ONLY))
chin_GTseq_ONLY <- chin_GTseq_ONLY[-1,]
colnames(chin_GTseq_ONLY) <- columns
chin_GTseq_ONLY <- cbind(SampleName, chin_GTseq_ONLY)
rownames(chin_GTseq_ONLY) <- c()
chin_GTseq_ONLY[1:6,1:6]

#a population map that has all the individuals and the populations that they belong to. 
#this file has been edited to combined the Keek pops 
chin_GTseq_ONLY_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/PopulationMap_GTseqDataONLY.txt", header =TRUE)
head(chin_GTseq_ONLY_popINFO)
dim(chin_GTseq_ONLY_popINFO)

#Duplicated loci
GTseqduplicates <- read.table("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/DuplicateLoci/polyGenResults_haplotypes_visualDuplicates_genoCode.txt", header= TRUE)
GTseqduplicates[1:5,1:10]

################### GTSEQ only##############
chin_GTseq_ONLY_rmvdup <- chin_GTseq_ONLY
chin_GTseq_ONLY_rmvdup[1:5,1:5]
dim(chin_GTseq_ONLY_rmvdup)

###remove known duplicated loci from the genepop file
duplicate_loci <- rownames(GTseqduplicates)
head(duplicate_loci)
length(duplicate_loci)

chin_GTseq_ONLY_rmvdup <- chin_GTseq_ONLY_rmvdup[,!(colnames(chin_GTseq_ONLY_rmvdup)%in% duplicate_loci)]
dim(chin_GTseq_ONLY_rmvdup)
chin_GTseq_ONLY_rmvdup[1:5,1:15]
GTseq_ONLYnodups_nloci <- colnames(chin_GTseq_ONLY_rmvdup)[-1]

####figure out which of the genotypes in the  GTseq only data set are more than 4 characters long ###
test_set <- chin_GTseq_ONLY_rmvdup[6,-1] #tested a couple different rows

list_to_dump <- vector("list",length(test_set))

for (i in 1:length(test_set)){
  test_a <- as.character(test_set[[i]])
  if (nchar(test_a) >4) {
    list_to_dump[[i]] <- colnames(test_set[i])
  }else{
    next()
  }
}

list_to_dump <- list_to_dump[-which(sapply(list_to_dump, is.null))] #remove the NULLS
length(list_to_dump)

#### Remove the loci that were found to be duplicated 
chin_GTseq_ONLY_rmvdup2 <- chin_GTseq_ONLY_rmvdup[,!(colnames(chin_GTseq_ONLY_rmvdup)%in% list_to_dump)]
dim(chin_GTseq_ONLY_rmvdup2)
chin_GTseq_ONLY_rmvdup2[1:5,1:15]
chin_GTseq_ONLY_rmvdup2[,1021]<- NULL

###Remove the 847 loci
#change the format of the loci names, replace any - with . to match the format of chin_comb_RADGTseq_loci
mid_filter <- colnames(chin_GTseq_ONLY_rmvdup2)
mid_filter<-gsub("-",".",mid_filter)
colnames(chin_GTseq_ONLY_rmvdup2) <- mid_filter

chin_GTseq_ONLY_rmvdup2_rmv847 <- chin_GTseq_ONLY_rmvdup2[,!(colnames(chin_GTseq_ONLY_rmvdup2)%in% chin_comb_RADGTseq_loci)]
dim(chin_GTseq_ONLY_rmvdup2_rmv847)
chin_GTseq_ONLY_rmvdup2_rmv847[1:5,1:15]
colnames(chin_GTseq_ONLY_rmvdup2_rmv847)

# 
# outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/chin_GTseq_ONLY_rmvdup2_rmv847", "wb")
# write.table(chin_GTseq_ONLY_rmvdup2_rmv847,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
# close(outputFile)

############# figure out which other loci are not able to have MAF calculated#
dim(chin_GTseq_ONLY_rmvdup2_rmv847)
chin_GTseq_ONLY_rmvdup2_rmv847[1:5,1:15]

list_to_look_at <- vector("list",180)

for (i in 2:dim(chin_GTseq_ONLY_rmvdup2_rmv847)[2]){
  templocus<-as.character(chin_GTseq_ONLY_rmvdup2_rmv847[,i])
  countalleles<-unique(templocus)
  if (length(countalleles) > 5){
    list_to_look_at[[i]] <- colnames(chin_GTseq_ONLY_rmvdup2_rmv847[i])
  } else{
    next()
  }
}

list_to_look_at <- list_to_look_at[-which(sapply(list_to_look_at, is.null))] #remove the NULLS
length(list_to_look_at)

chin_GTseq_ONLY_rmvdup2_rmv847_2 <- chin_GTseq_ONLY_rmvdup2_rmv847[,!(colnames(chin_GTseq_ONLY_rmvdup2_rmv847)%in% list_to_look_at)]
dim(chin_GTseq_ONLY_rmvdup2_rmv847_2)
chin_GTseq_ONLY_rmvdup2_rmv847_2[1:5,1:15]
colnames(chin_GTseq_ONLY_rmvdup2_rmv847_2)

###### Merge the Pop info with the genotypes ######
#assign each sample to a population 
chin_GTseq_pop<-merge(chin_GTseq_ONLY_popINFO,chin_GTseq_ONLY_rmvdup2_rmv847_2, by="SampleName")
chin_GTseq_pop[1:5,1:10]
dim(chin_GTseq_pop)

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


######### Test MAF >= 0.05 in at least one of the 4 main regions: BroadReportingGroup (BroadRegion, column 5) ############
colnames(chin_GTseq_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_GTseq_pop$Broad_Region))
nloci <- length(colnames(chin_GTseq_ONLY_rmvdup2_rmv847_2))-1

GTseq_BREG_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(GTseq_BREG_MAF)<-as.vector(unique(chin_GTseq_pop$Broad_Region)) #name the rows by the population names
colnames(GTseq_BREG_MAF)<-colnames(chin_GTseq_ONLY_rmvdup2_rmv847_2)[-1]
dim(GTseq_BREG_MAF)
GTseq_BREG_MAF[1:3,1:6]

##### Get the MAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_GTseq_pop$Broad_Region)[i])
  tempGeno<-subset(chin_GTseq_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(7:dim(tempGeno)[2])],2,calculateMAF)
  GTseq_BREG_MAF[i,]<-c(tempMAF)
}

unique(chin_GTseq_pop[,18])
any(is.na(chin_GTseq_pop))

GTseq_BREG_MAF[1:3,1:5]
dim(GTseq_BREG_MAF)

#colnames(GTseq_BREG_MAF[,1])





# #######################################Havent gotten this far yet


# outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/chin_GTseq_pop.txt", "wb")
# write.table(chin_GTseq_pop,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_BREGMAF<-vector()
loci_BREGMAF<-apply(GTseq_BREG_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
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















