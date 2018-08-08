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
chin_comb_RADGTseq <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", sep="", header = TRUE, colClasses="factor", check.names=FALSE)
dim(chin_comb_RADGTseq)
chin_comb_RADGTseq[1:6,1:6]

##Import the list of the 847 loci names::
RADGTseq_loci <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/List_of_847_loci_names.txt", header = TRUE)
dim(RADGTseq_loci)
head(RADGTseq_loci)

#load the mixed one snp and haplotype genepop file that is all the 1092 GTSEQ ONLY data
chin_GTseq_ONLY <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/combinedPrimerProbe_amplicon_polyResults_genepop.txt", sep="", header = TRUE, colClasses="factor", row.names = 1)
dim(chin_GTseq_ONLY)
chin_GTseq_ONLY[1:6,1:6]

#Transpose the data so that each individual is a row
chin_GTseq_ONLY <- as.data.frame(t(chin_GTseq_ONLY))
dim(chin_GTseq_ONLY)
chin_GTseq_ONLY[1:10,1080:1092]

SampleName <-rownames(chin_GTseq_ONLY)
length(SampleName)
chin_GTseq_ONLY <- cbind(SampleName, chin_GTseq_ONLY)
colnames(chin_GTseq_ONLY)
rownames(chin_GTseq_ONLY) <- c()
chin_GTseq_ONLY[1:6,1:6]
chin_GTseq_ONLY[1:10,1080:1093]
dim(chin_GTseq_ONLY)

# ##########write this transposed file out to use as the genepop file
# outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/combinedPrimerProbe_amplicon_polyResults_genepop_transposed.txt", "wb")
# write.table(chin_GTseq_ONLY,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n", sep= "\t")
# close(outputFile)

#a population map that has all the individuals and the populations that they belong to. 
#this file has been edited to combined the Keek pops 
chin_GTseq_ONLY_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/PopulationMap_GTseqDataONLY.txt", header =TRUE)
head(chin_GTseq_ONLY_popINFO)
dim(chin_GTseq_ONLY_popINFO)

#Duplicated loci
GTseqduplicates_raw <- read.table("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/DuplicateLoci/polyGenResults_haplotypes_visualDuplicates_genoCode.txt", colClasses="factor", header= TRUE)
GTseqduplicates_raw[1:5,1:10]
dim(GTseqduplicates_raw)

## List of the haplotypes. I used the allelekey- produced from Garrett's perl script to take a haplotype file in to a Genepop file
#for identifying haplotypes
allele_key <-read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/convertHaplotypes_GTseqonly_hapLoci.txt", header =FALSE)
colnames(allele_key) <- "Haplotypes"
head(allele_key)
dim(allele_key)

#######################################################################################################
################### GTSEQ only#########################################################################


#######################################################################################################
###################Sort the data into duplicated and non duplicated ###################################
chin_GTseq_ONLY_rmvdup <- chin_GTseq_ONLY
chin_GTseq_ONLY_rmvdup[1:5,1:5]
dim(chin_GTseq_ONLY_rmvdup)

###remove known duplicated loci from the non duplicated data set 
duplicate_loci <- rownames(GTseqduplicates_raw)
head(duplicate_loci)
length(duplicate_loci)

chin_GTseq_ONLY_rmvdup <- chin_GTseq_ONLY_rmvdup[,!(colnames(chin_GTseq_ONLY_rmvdup)%in% duplicate_loci)]
dim(chin_GTseq_ONLY_rmvdup)
chin_GTseq_ONLY_rmvdup[1:5,1:15]
GTseq_ONLYnodups_nloci <- colnames(chin_GTseq_ONLY_rmvdup)[-1]

####figure out which of the genotypes in the  GTseq only data set are more than 4 characters long ###
##to remove them from the duplicated set and add them to the duplicated genotypes
test_set <- chin_GTseq_ONLY_rmvdup[25,-1] #tested a couple different rows
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
list_to_dump

#verify the genotypes of the loci in the list. 
# chin_GTseq_ONLY[,"RAD4642"]
# "Ots_126619-400"
# "Ots_U2305-63"
# "RAD32242"
# "RAD34754"
# "RAD44889"
# "RAD4642"


####create a data set that is the duplicated loci  ##################
dim(GTseqduplicates_raw) 

#Transpose the data so that each individual is a row
GTseqduplicates <- as.data.frame(t(GTseqduplicates_raw))
dim(GTseqduplicates)
GTseqduplicates[1:10,1:10]

SampleName <-rownames(GTseqduplicates)
length(SampleName)
GTseqduplicates <- cbind(SampleName, GTseqduplicates)
rownames(GTseqduplicates) <- c()
dim(GTseqduplicates)

###include the loci that were found to be duplicated that werent in the duplicated list######
chin_GTseq_ONLY_rmvdup_list_to_dump <- chin_GTseq_ONLY[,((colnames(chin_GTseq_ONLY)%in% list_to_dump) | (colnames(chin_GTseq_ONLY) == "SampleName"))]
dim(chin_GTseq_ONLY_rmvdup_list_to_dump)
chin_GTseq_ONLY_rmvdup_list_to_dump[1:5,1:6]

chin_GTseq_ONLY_duplicated<- merge(chin_GTseq_ONLY_rmvdup_list_to_dump, GTseqduplicates, by = "SampleName")
dim(chin_GTseq_ONLY_duplicated)
chin_GTseq_ONLY_duplicated[1:5,1:6]

####verify that all the genotypes in the duplicated set are 8 characters long  ###
test_set_dup <- chin_GTseq_ONLY_duplicated[25,-1] #tested a couple different rows
list_to_dup_dump <- vector("list",length(test_set_dup))

for (i in 1:length(test_set_dup)){
  test_a <- as.character(test_set_dup[[i]])
  if (nchar(test_a) <8) {
    list_to_dup_dump[[i]] <- colnames(test_set_dup[i])
  }else{
    next()
  }
}

list_to_dup_dump <- list_to_dup_dump[-which(sapply(list_to_dup_dump, is.null))] #remove the NULLS
length(list_to_dup_dump)
list_to_dup_dump

#### Remove the loci that were found to be duplicated from the nonduplicated dataset 
chin_GTseq_ONLY_rmvdup2 <- chin_GTseq_ONLY_rmvdup[,!(colnames(chin_GTseq_ONLY_rmvdup)%in% list_to_dump)]
dim(chin_GTseq_ONLY_rmvdup2)
chin_GTseq_ONLY_rmvdup2[1:5,1:15]

####Identify and remove loci that are just 0000 from the non duplicated set
testchin_GTseq_ONLY_rmvdup2<- chin_GTseq_ONLY_rmvdup2[,2:(dim(chin_GTseq_ONLY_rmvdup2)[2])]
dim(testchin_GTseq_ONLY_rmvdup2)
testchin_GTseq_ONLY_rmvdup2[1:5,1:15]

list_to_dump_0000 <- vector("list",length(testchin_GTseq_ONLY_rmvdup2))
#testchin_GTseq_ONLY_rmvdup2[,"RAD83732"]

for (i in 1:dim(testchin_GTseq_ONLY_rmvdup2)[2]){
  test_column <- testchin_GTseq_ONLY_rmvdup2[,i]
  unique_genotypes <- as.vector(unique(test_column))
  if ((length(unique_genotypes) == 1) & (unique_genotypes == "0000" )) { #this throws a warning, but it's OK, if the first condition is met, there is only one option for the second condition
    list_to_dump_0000[[i]] <- colnames(testchin_GTseq_ONLY_rmvdup2[i])
  }else{
    #print(paste(colnames(testchin_GTseq_ONLY_rmvdup2[i])))
    next()
  }
}

list_to_dump_0000 <- list_to_dump_0000[-which(sapply(list_to_dump_0000, is.null))] #remove the NULLS
length(list_to_dump_0000)
list_to_dump_0000

#### Remove the loci that were  0000
chin_GTseq_ONLY_rmvdup2 <- chin_GTseq_ONLY_rmvdup2[,!(colnames(chin_GTseq_ONLY_rmvdup2)%in% list_to_dump_0000)]
dim(chin_GTseq_ONLY_rmvdup2)
chin_GTseq_ONLY_rmvdup2[1:5,1:15]

### Remove the loci that are also RADGTseq_loci from the non duplicated set
sum(colnames(chin_GTseq_ONLY_rmvdup2)%in% RADGTseq_loci$X847_Loci)
chin_GTseq_ONLY_rmvdup2_rmv847 <- chin_GTseq_ONLY_rmvdup2[,!(colnames(chin_GTseq_ONLY_rmvdup2)%in% RADGTseq_loci$X847_Loci)]
dim(chin_GTseq_ONLY_rmvdup2_rmv847)
chin_GTseq_ONLY_rmvdup2_rmv847[1:5,1:15]
colnames(chin_GTseq_ONLY_rmvdup2_rmv847)

#Remove the loci that are also in the RADGTseq_loci from the duplicated set (verify that there arent any)
sum(colnames(chin_GTseq_ONLY_duplicated)%in% RADGTseq_loci$X847_Loci)



#######################################################################################################
#### Split nonduplicated genotypes into Haplotype only data set. ######################################

chin_comb_GTseq_Hap <- chin_GTseq_ONLY_rmvdup2_rmv847[,colnames(chin_GTseq_ONLY_rmvdup2_rmv847) %in% allele_key$Haplotypes]
dim(chin_comb_GTseq_Hap)
SampleName <- chin_GTseq_ONLY_rmvdup2_rmv847[,1]
chin_comb_GTseq_Hap <- cbind(SampleName, chin_comb_GTseq_Hap)#merge the sample names back on to the genotypes
chin_comb_GTseq_Hap[1:5,1:5]

##### Split non duplicated genotypes into OneSNP only data set. 
chin_comb_GTseq_OneSNP <- chin_GTseq_ONLY_rmvdup2_rmv847[,!(colnames(chin_GTseq_ONLY_rmvdup2_rmv847) %in% allele_key$Haplotypes)]
dim(chin_comb_GTseq_OneSNP)
chin_comb_GTseq_OneSNP[1:5,1:5]

####### Merge the Pop info with the non duplicated genotypes ######
#assign each sample to a population 
#haplotype data
chin_comb_GTseq_Hap_pop<-merge(chin_GTseq_ONLY_popINFO,chin_comb_GTseq_Hap, by="SampleName")
chin_comb_GTseq_Hap_pop[1:5,1:15]
dim(chin_comb_GTseq_Hap_pop)

#ONE snp per tag data
chin_comb_GTseq_OneSNP_pop<-merge(chin_GTseq_ONLY_popINFO,chin_comb_GTseq_OneSNP, by="SampleName")
chin_comb_GTseq_OneSNP_pop[1:5,1:15]
dim(chin_comb_GTseq_OneSNP_pop)


#######################################################################################################
##############Minor Allele Frequency MinAF#############################################################
#ONLY For the one SNP per tag data set

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

######### Test MinAF >= 0.05 in at least one of the 3 main regions: BroadReportingGroup (Broad_Region, column 5) ############
colnames(chin_comb_GTseq_OneSNP_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_comb_GTseq_OneSNP_pop$Broad_Region))
nloci <- length(colnames(chin_comb_GTseq_OneSNP_pop))-6

GTseq_BREG_MinAF<-matrix(nrow=npops, ncol=nloci) 
rownames(GTseq_BREG_MinAF)<-as.vector(unique(chin_comb_GTseq_OneSNP_pop$Broad_Region)) #name the rows by the population names
colnames(GTseq_BREG_MinAF)<-colnames(chin_comb_GTseq_OneSNP_pop)[-c(1:6)]
dim(GTseq_BREG_MinAF)
GTseq_BREG_MinAF[1:3,1:6]

##### Get the MinAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_comb_GTseq_OneSNP_pop$Broad_Region)[i])
  tempGeno<-subset(chin_comb_GTseq_OneSNP_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(7:dim(tempGeno)[2])],2,calculateMAF)
  GTseq_BREG_MinAF[i,]<-c(tempMAF)
}

GTseq_BREG_MinAF[1:3,1:5]
dim(GTseq_BREG_MinAF)

#test for specific loci
#chin_comb_GTseq_OneSNP_pop[,"Ots_106313-729"]

##### Filter the loci by the MinAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_BREGMinAF<-vector()
loci_BREGMinAF<-apply(GTseq_BREG_MinAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
loci_BREGMinAF<-data.frame(keynames=names(loci_BREGMinAF), value=loci_BREGMinAF, row.names = NULL)
colnames(loci_BREGMinAF)<-c("Locus", "BroadRegionMAF")
head(loci_BREGMinAF)

#if the loci does not pass the test, delete it
locipassed_BREGMinAF<-vector()
locipassed_BREGMinAF<-subset(loci_BREGMinAF, loci_BREGMinAF[,2]!=0)
dim(locipassed_BREGMinAF)
head(locipassed_BREGMinAF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMinAF <- ((dim(locipassed_BREGMinAF)[1])/149)*100
Percent_locipassed_BREGMinAF

#get list of loci that made it past this filter
locipassed_BREGMinAF_keepList <- locipassed_BREGMinAF$Locus
length(locipassed_BREGMinAF_keepList)
head(locipassed_BREGMinAF_keepList)
drop.levels(locipassed_BREGMinAF_keepList)

###output a list of the loci that passed the minor allele filter for the one snp no duplicated set
# outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/GT_ONESNP_locipassed_BREGMinAF_keepList.txt", "wb")
# write.table(locipassed_BREGMinAF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

##make a genotype dataset for the loci that passed the one snp non duplicated minor allele filter
chin_comb_GTseq_OneSNP_keepGenos <- chin_GTseq_ONLY_rmvdup2_rmv847[,(colnames(chin_GTseq_ONLY_rmvdup2_rmv847) %in% locipassed_BREGMinAF_keepList)| (colnames(chin_GTseq_ONLY_rmvdup2_rmv847) == "SampleName")]
dim(chin_comb_GTseq_OneSNP_keepGenos)
chin_comb_GTseq_OneSNP_keepGenos[1:5,1:5]

###output a genotypes table for the loci that passed the minor allele filter for the one snp no duplicated set
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/chin_comb_GTseq_OneSNP_keepGenos.txt", "wb")
write.table(chin_comb_GTseq_OneSNP_keepGenos,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)


#######################################################################################################
############## Major Allele Frequency MajAF############################################################

#ONLY For the Haplotype data set
chin_comb_GTseq_Hap_pop[1:10,1:15]
dim(chin_comb_GTseq_Hap_pop)
dim(chin_comb_GTseq_Hap_pop)[2]-6


#Function to get major allele frequency for each locus #this can work with Haplotypes with up to 7 alleles, look at how many you have in the ALLELEKey.txt
calculateMAJF<-function(genotypes){
  genotypeList<-sort(unique(genotypes))
  genotypeList<-genotypeList[genotypeList != "0000"]
  allelesList1<-substr(genotypeList,1,2)
  allelesList2<-substr(genotypeList,3,4)
  allelesList<-unique(c(allelesList1,allelesList2))
  allele1Counts<-sum(str_count(genotypes,allelesList[1]))
  allele2Counts<-sum(str_count(genotypes,allelesList[2]))
  allele3Counts<-sum(str_count(genotypes,allelesList[3]))
  allele4Counts<-sum(str_count(genotypes,allelesList[4]))
  allele5Counts<-sum(str_count(genotypes,allelesList[5]))
  allele6Counts<-sum(str_count(genotypes,allelesList[6]))
  allele7Counts<-sum(str_count(genotypes,allelesList[7]))
  options <- c(allele1Counts, allele2Counts ,allele3Counts ,allele4Counts ,allele5Counts ,allele6Counts ,allele7Counts)
  max_options <- max(options, na.rm= TRUE)
  total <- sum(options, na.rm = TRUE)
  MAJF<-max_options/total
  return(MAJF)
}


######### Test MajAF >= 0.05 in at least one of the 3 main regions: Broad_Region (column 5) ############
colnames(chin_comb_GTseq_Hap_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_comb_GTseq_Hap_pop$Broad_Region))
nloci <- length(colnames(chin_comb_GTseq_Hap_pop))-6

GTseq_BREG_MajAF<-matrix(nrow=npops, ncol=nloci) 
rownames(GTseq_BREG_MajAF)<-as.vector(unique(chin_comb_GTseq_Hap_pop$Broad_Region)) #name the rows by the population names
colnames(GTseq_BREG_MajAF)<-colnames(chin_comb_GTseq_Hap_pop[-c(1:6)])
dim(GTseq_BREG_MajAF)
GTseq_BREG_MajAF[1:3,1:6]

##### Get the MAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_comb_GTseq_Hap_pop$Broad_Region)[i])
  tempGeno<-subset(chin_comb_GTseq_Hap_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(7:dim(tempGeno)[2])],2,calculateMAJF)
  GTseq_BREG_MajAF[i,]<-c(tempMAF)
}

GTseq_BREG_MajAF[1:3,1:5]
dim(GTseq_BREG_MajAF)


##### Filter the loci by the MajAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_BREGMAJF<-vector()
loci_BREGMAJF<-apply(GTseq_BREG_MajAF,2,function(x) sum(x <=0.95, na.rm=TRUE))
loci_BREGMAJF<-data.frame(keynames=colnames(GTseq_BREG_MajAF), value=loci_BREGMAJF, row.names = NULL)
colnames(loci_BREGMAJF)<-c("Locus", "BroadRegionMAJF")
head(loci_BREGMAJF)

#if the loci does not pass the test, delete it
locipassed_BREGMAJF<-vector()
locipassed_BREGMAJF<-subset(loci_BREGMAJF, loci_BREGMAJF[,2]!=0)
dim(locipassed_BREGMAJF)
head(locipassed_BREGMAJF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMAJF <- ((dim(locipassed_BREGMAJF)[1])/26)*100
Percent_locipassed_BREGMAJF

#get list of loci that made it past this filter
locipassed_BREGMAJF_keepList <- locipassed_BREGMAJF$Locus
length(locipassed_BREGMAJF_keepList)
head(locipassed_BREGMAJF_keepList)
locipassed_BREGMAJF_keepList <- drop.levels(locipassed_BREGMAJF_keepList)


################## Output lists of Loci that passed MajAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/GTseq_Haplocipassed_BREGMAJF_keepList.txt", "wb")
write.table(locipassed_BREGMAJF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

#make a genotype data set for the loci that passed the major allele frequency filter in the haplotype non duplicated data set
chin_comb_GTseq_Hap_keepGenos <- chin_comb_GTseq_Hap_pop[,(colnames(chin_comb_GTseq_Hap_pop) %in% locipassed_BREGMAJF_keepList)| (colnames(chin_comb_GTseq_Hap_pop) == "SampleName")]
dim(chin_comb_GTseq_Hap_keepGenos)
chin_comb_GTseq_Hap_keepGenos[1:5,1:5]

################## Output the genotypes for the non duplicated loci that passed MajAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/chin_comb_GTseq_Hap_keepGenos.txt", "wb")
write.table(chin_comb_GTseq_Hap_keepGenos,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

# make a dataset combines all the loci that passed the two filters,  for the non duplicated loci, shuld all have 4 character coding###################
GTseqONLY_nonduplicated_genos_keepLoci <- merge (chin_comb_GTseq_OneSNP_keepGenos, chin_comb_GTseq_Hap_keepGenos, by = "SampleName")
dim(GTseqONLY_nonduplicated_genos_keepLoci)
GTseqONLY_nonduplicated_genos_keepLoci[1:5,1:7]

################## Output the combined genotypes for the non duplicated loci that passed either filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/GTseqONLY_nonduplicated_genos_keep.txt", "wb")
write.table(GTseqONLY_nonduplicated_genos_keepLoci,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

#######################################################################################################
############################Look at the allele frequencies of the duplicated loci######################

#This is the full duplicated genotypes data set
dim(chin_GTseq_ONLY_duplicated)

####### Merge the Pop info with the genotypes ######
chin_GTseq_ONLY_duplicated_pop<-merge(chin_GTseq_ONLY_popINFO,chin_GTseq_ONLY_duplicated, by="SampleName")
chin_GTseq_ONLY_duplicated_pop[1:5,1:15]
dim(chin_GTseq_ONLY_duplicated_pop)

######### Calculate the allele Frequency for each of the duplicated loci

calculateAlleleFreq<-function(genotypes){ #this gives all of them in a list 
  locus_inQ <- genotypes
  genotypeList<-locus_inQ[locus_inQ != "00000000"]
  tabled_set <- sort(table(genotypeList))
  genotypes<-names(tabled_set)
  allelesList1<-substr(genotypeList,1,2)
  allelesList2<-substr(genotypeList,3,4)
  allelesList3<-substr(genotypeList,5,6)
  allelesList4<-substr(genotypeList,7,8)
  alleles<-unique(c(allelesList1,allelesList2,allelesList3,allelesList4))
  alleleCounts<-unlist(lapply(alleles,function(x) sum(str_count(genotypeList,x))))
  alleleFreqs<-alleleCounts/sum(alleleCounts)
  names(alleleFreqs)<-alleles
  return(alleleFreqs)
}

##### Get the Allele frequency for each locus by region 
GTseq_duplicate_AF <- list()

for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_GTseq_ONLY_duplicated_pop$Broad_Region)[i])
  tempGeno<-subset(chin_GTseq_ONLY_duplicated_pop, Broad_Region == tempRegion)
  tempAF<-apply(tempGeno[,c(7:dim(tempGeno)[2])],2,calculateAlleleFreq)
  GTseq_duplicate_AF[[i]]<-c(tempAF)
}

GTseq_duplicate_AF
length(GTseq_duplicate_AF)

######### Calculate the major allele Frequency for each of the duplicated loci
#this can work with Haplotypes with up to 7 alleles, look at how many you have in the ALLELEKey.txt

calculateMAJF_dup<-function(genotypes){
  genotypeList<-sort(unique(genotypes))
  genotypeList<-genotypeList[genotypeList != "00000000"]
  allelesList1<-substr(genotypeList,1,2)
  allelesList2<-substr(genotypeList,3,4)
  allelesList3<-substr(genotypeList,5,6)
  allelesList4<-substr(genotypeList,7,8)
  allelesList<-unique(c(allelesList1,allelesList2,allelesList3, allelesList4))
  allele1Counts<-sum(str_count(genotypes,allelesList[1]))
  allele2Counts<-sum(str_count(genotypes,allelesList[2]))
  allele3Counts<-sum(str_count(genotypes,allelesList[3]))
  allele4Counts<-sum(str_count(genotypes,allelesList[4]))
  allele5Counts<-sum(str_count(genotypes,allelesList[5]))
  allele6Counts<-sum(str_count(genotypes,allelesList[6]))
  allele7Counts<-sum(str_count(genotypes,allelesList[7]))
  options <- c(allele1Counts, allele2Counts ,allele3Counts ,allele4Counts ,allele5Counts ,allele6Counts ,allele7Counts)
  max_options <- max(options, na.rm= TRUE)
  total <- sum(options, na.rm = TRUE)
  MAJF_dup<-max_options/total
  return(MAJF_dup)
}

##### Get the major allele frequency for each locus by region 
#create a dataframe for our MAF results
npops<- length(unique(chin_GTseq_ONLY_duplicated_pop$Broad_Region))
nloci <- length(colnames(chin_GTseq_ONLY_duplicated_pop))-6

GTseq_duplicate_MAJF<-matrix(nrow=npops, ncol=nloci) 
rownames(GTseq_duplicate_MAJF)<-as.vector(unique(chin_GTseq_ONLY_duplicated_pop$Broad_Region)) #name the rows by the population names
colnames(GTseq_duplicate_MAJF)<-colnames(chin_GTseq_ONLY_duplicated_pop[-c(1:6)])
dim(GTseq_duplicate_MAJF)
GTseq_duplicate_MAJF[1:3,1:6]

for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_GTseq_ONLY_duplicated_pop$Broad_Region)[i])
  tempGeno<-subset(chin_GTseq_ONLY_duplicated_pop, Broad_Region == tempRegion)
  tempAF<-apply(tempGeno[,c(7:dim(tempGeno)[2])],2,calculateMAJF_dup)
  GTseq_duplicate_MAJF[i,]<-c(tempAF)
}

GTseq_duplicate_MAJF

##### Filter the loci by the MinAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_dup_BREGMAJF<-vector()
loci_dup_BREGMAJF<-apply(GTseq_duplicate_MAJF,2,function(x) sum(x <=0.95, na.rm=TRUE))
loci_dup_BREGMAJF<-data.frame(keynames=colnames(GTseq_duplicate_MAJF), value=loci_dup_BREGMAJF, row.names = NULL)
colnames(loci_dup_BREGMAJF)<-c("Locus", "BroadRegionMAJF")
head(loci_dup_BREGMAJF)
dim(loci_dup_BREGMAJF)

#if the loci does not pass the test, delete it
locipassed_BREGMAJF<-vector()
locipassed_BREGMAJF<-subset(loci_dup_BREGMAJF, loci_dup_BREGMAJF[,2]!=0)
dim(locipassed_BREGMAJF)
head(locipassed_BREGMAJF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMAJF <- ((dim(locipassed_BREGMAJF)[1])/72)*100
Percent_locipassed_BREGMAJF

#get list of loci that made it past this filter
locipassed_BREGMAJF_keepList <- locipassed_BREGMAJF$Locus
length(locipassed_BREGMAJF_keepList)
head(locipassed_BREGMAJF_keepList)
locipassed_BREGMAJF_keepList <- drop.levels(locipassed_BREGMAJF_keepList)

#subset the genotypes of the loci that passed this filter

chin_GTseq_ONLY_duplicated_pop_BREGMAJF_keep <- chin_GTseq_ONLY_duplicated_pop[,((colnames(chin_GTseq_ONLY_duplicated_pop) %in% locipassed_BREGMAJF_keepList) | (colnames(chin_GTseq_ONLY_duplicated_pop) == "SampleName"))]
dim(chin_GTseq_ONLY_duplicated_pop_BREGMAJF_keep)
chin_GTseq_ONLY_duplicated_pop_BREGMAJF_keep[1:5,1:7]


################## Output the genotypes of the Loci that passed MajAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/chin_GTseq_ONLY_duplicated_BREGMAJF_genos_keep.txt", "wb")
write.table(chin_GTseq_ONLY_duplicated_pop_BREGMAJF_keep,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

