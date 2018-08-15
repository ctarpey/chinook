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

##There was a naming convention issue, so I remade the genepop file, this is the upated version
#load the mixed haplotype and non haplotype genepop file that is all the 847 Combined RADGTseq data
chin_comb_RADGTseq <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_Rnew.txt", sep="", header = TRUE, colClasses="factor", check.names=FALSE)
dim(chin_comb_RADGTseq)
chin_comb_RADGTseq[1:6,1:6]

# #load the mixed haplotype and non haplotype genepop file that is all the 847 Combined RADGTseq data- ORIGINAL VERSION
# chin_comb_RADGTseq <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADtaqAmp_Rcombined_genepop_R.txt", sep="", header = TRUE, colClasses="factor")
# dim(chin_comb_RADGTseq)
# chin_comb_RADGTseq[1:6,1:6]

#prep for later- the names of the loci and the number: 
chin_comb_RADGTseq_loci <- colnames(chin_comb_RADGTseq)[2:(dim(chin_comb_RADGTseq)[2])]
head(chin_comb_RADGTseq_loci)
nloci <- length(chin_comb_RADGTseq_loci)
nloci

#a population map that has all the individuals and the populations that they belong to. 
#this file has been edited to combined the Keek pops 
chin_RADGTseq_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/PopulationMap_847RAD_taq_genepop.txt", header =TRUE)
head(chin_RADGTseq_popINFO)
dim(chin_RADGTseq_popINFO)

## List of the haplotypes. I used the allelekey- produced from Garrett's perl script to take a haplotype file in to a Genepop file
#for identifying haplotypes

allele_key <-read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/convertHaplotypes_RADGTseq_hapLoci.txt", header =FALSE)
colnames(allele_key) <- "Haplotypes"
head(allele_key)
dim(allele_key)


###########################Combined RAD GTseq starting with 847 loci #####################

##### Split genotypes into  Haplotype only data set. 
chin_comb_RADGTseq_Hap <- chin_comb_RADGTseq[,colnames(chin_comb_RADGTseq) %in% allele_key$Haplotypes]
dim(chin_comb_RADGTseq_Hap)
SampleName <- chin_comb_RADGTseq[,1]
chin_comb_RADGTseq_Hap <- cbind(SampleName, chin_comb_RADGTseq_Hap)#merge the sample names back on to the genotypes
chin_comb_RADGTseq_Hap[1:5,1:5]

##### Split genotypes into OneSNP only data set. 
chin_comb_RADGTseq_OneSNP <- chin_comb_RADGTseq[,!(colnames(chin_comb_RADGTseq) %in% allele_key$Haplotypes)]
dim(chin_comb_RADGTseq_OneSNP)
chin_comb_RADGTseq_OneSNP[1:5,1:5]

####### Merge the Pop info with the genotypes ######
#assign each sample to a population 
#haplotype data
chin_comb_RADGTseq_Hap_pop<-merge(chin_RADGTseq_popINFO,chin_comb_RADGTseq_Hap, by="SampleName")
chin_comb_RADGTseq_Hap_pop[1:5,1:15]
dim(chin_comb_RADGTseq_Hap_pop)

#ONE snp per tag data
chin_comb_RADGTseq_OneSNP_pop<-merge(chin_RADGTseq_popINFO,chin_comb_RADGTseq_OneSNP, by="SampleName")
chin_comb_RADGTseq_OneSNP_pop[1:5,1:15]
dim(chin_comb_RADGTseq_OneSNP_pop)

##############Minor Allele Frequency MinAF#############
#ONLY For the one SNP per tag data set

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

chin_comb_RADGTseq_OneSNP_pop[1:10,1:15]
dim(chin_comb_RADGTseq_OneSNP_pop)[2]-5

######### Test MinAF >= 0.05 in at least one of 17 pops ############
#create a dataframe for our MAF results
npops<- length(unique(chin_comb_RADGTseq_OneSNP_pop$Population))
nloci <- dim(chin_comb_RADGTseq_OneSNP_pop)[2]-5 # the 5 is the number of columns of population data on the front
loci <- colnames(chin_comb_RADGTseq_OneSNP_pop)[c(6:dim(chin_comb_RADGTseq_OneSNP_pop)[2])]
length(loci)

RADGTseq_POP_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(RADGTseq_POP_MAF)<-as.vector(unique(chin_comb_RADGTseq_OneSNP_pop$Population)) #name the rows by the population names
colnames(RADGTseq_POP_MAF)<-loci
RADGTseq_POP_MAF[1:6,1:6]
dim(RADGTseq_POP_MAF)

#### Get the MinAF for each locus by population 
for (i in 1:17){ #17 is the unique number of populations we have
  tempPop<-(unique(chin_comb_RADGTseq_OneSNP_pop$Population)[i])
  tempGeno<-subset(chin_comb_RADGTseq_OneSNP_pop, Population == tempPop)
  tempGeno[,1:10]
  # dim(tempGeno)
  tempMAF<-apply(tempGeno[,c(6:dim(chin_comb_RADGTseq_OneSNP_pop)[2])],2,calculateMAF)
  RADGTseq_POP_MAF[i,]<-c(tempMAF)
}

RADGTseq_POP_MAF[1:17,1:5]
dim(RADGTseq_POP_MAF)
#colnames(RADGTseq_POP_MAF[,1])

##### Filter the loci by the MinAF results by population
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

######### Test MinAF >= 0.05 in at least one of the 3 main regions: BroadReportingGroup (BroadRegion, column 5) ############
colnames(chin_comb_RADGTseq_OneSNP_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_comb_RADGTseq_OneSNP_pop$Broad_Region))

RADGTseq_BREG_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(RADGTseq_BREG_MAF)<-as.vector(unique(chin_comb_RADGTseq_OneSNP_pop$Broad_Region)) #name the rows by the population names
colnames(RADGTseq_BREG_MAF)<-loci
dim(RADGTseq_BREG_MAF)
RADGTseq_BREG_MAF[1:3,1:6]

##### Get the MinAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_comb_RADGTseq_OneSNP_pop$Broad_Region)[i])
  tempGeno<-subset(chin_comb_RADGTseq_OneSNP_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(6:dim(chin_comb_RADGTseq_OneSNP_pop)[2])],2,calculateMAF)
  RADGTseq_BREG_MAF[i,]<-c(tempMAF)
}

RADGTseq_BREG_MAF[1:3,1:5]
dim(RADGTseq_BREG_MAF)
#colnames(RADGTseq_BREG_MAF[,1])

##### Filter the loci by the MinAF results by population
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
Percent_locipassed_BREGMAF <- ((dim(locipassed_BREGMAF)[1])/702)*100
Percent_locipassed_BREGMAF

#get list of loci that made it past this filter
locipassed_BREGMAF_keepList <- locipassed_BREGMAF$Locus
length(locipassed_BREGMAF_keepList)
head(locipassed_BREGMAF_keepList)
locipassed_BREGMAF_keepList <- drop.levels(locipassed_BREGMAF_keepList)


################# Filter out Loci that do not meet the MinAF threshold of 0.05 in 1 or more regions ##########################
filtered_MAF_Genos<-chin_comb_RADGTseq_OneSNP_pop[,((colnames(chin_comb_RADGTseq_OneSNP_pop) %in% locipassed_BREGMAF_keepList) | (colnames(chin_comb_RADGTseq_OneSNP_pop) %in% colnames(chin_RADGTseq_popINFO)))]
dim(filtered_MAF_Genos)
filtered_MAF_Genos[1:4,1:15]

################## Output lists of Loci that passed MinAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/locipassed_BroadREGMinAF_keepList.txt", "wb")
write.table(locipassed_BREGMAF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

chin_comb_RADGTseq_Hap_pop[1:5,1:5]

############## Major Allele Frequency MajAF#############
#ONLY For the Haplotype data set
chin_comb_RADGTseq_Hap_pop[1:10,1:15]
dim(chin_comb_RADGTseq_Hap_pop)[2]-5

######### Test MajAF >= 0.05 in at least one of the 3 main regions: Broad_Region (column 5) ############

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
colnames(chin_comb_RADGTseq_Hap_pop[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(chin_comb_RADGTseq_Hap_pop$Broad_Region))
nloci <- length(colnames(chin_comb_RADGTseq_Hap_pop))-5

RADGTseq_Hap_MAJF<-matrix(nrow=npops, ncol=nloci) 
rownames(RADGTseq_Hap_MAJF)<-as.vector(unique(chin_comb_RADGTseq_Hap_pop$Broad_Region)) #name the rows by the population names
colnames(RADGTseq_Hap_MAJF)<-colnames(chin_comb_RADGTseq_Hap_pop[-c(1:5)])
dim(RADGTseq_Hap_MAJF)
RADGTseq_Hap_MAJF[1:3,1:6]

##### Get the MAF for each locus by region 
for (i in 1:3){ #3 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(chin_comb_RADGTseq_Hap_pop$Broad_Region)[i])
  tempGeno<-subset(chin_comb_RADGTseq_Hap_pop, Broad_Region == tempRegion)
  tempMAF<-apply(tempGeno[,c(6:dim(tempGeno)[2])],2,calculateMAJF)
  RADGTseq_Hap_MAJF[i,]<-c(tempMAF)
}

RADGTseq_Hap_MAJF[1:3,1:5]
dim(RADGTseq_Hap_MAJF)

##### Filter the loci by the MajAF results by population
##retain loci that were at least 0.05 in any of the 4 broad regions
loci_BREGMAJF<-vector()
loci_BREGMAJF<-apply(RADGTseq_Hap_MAJF,2,function(x) sum(x <=0.95, na.rm=TRUE))
loci_BREGMAJF<-data.frame(keynames=colnames(RADGTseq_Hap_MAJF), value=loci_BREGMAJF, row.names = NULL)
colnames(loci_BREGMAJF)<-c("Locus", "BroadRegionMAJF")
head(loci_BREGMAJF)

#if the loci does not pass the test, delete it
locipassed_BREGMAJF<-vector()
locipassed_BREGMAJF<-subset(loci_BREGMAJF, loci_BREGMAJF[,2]!=0)
dim(locipassed_BREGMAJF)
head(locipassed_BREGMAJF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMAJF <- ((dim(locipassed_BREGMAJF)[1])/145)*100
Percent_locipassed_BREGMAJF

#get list of loci that made it past this filter
locipassed_BREGMAJF_keepList <- locipassed_BREGMAJF$Locus
length(locipassed_BREGMAJF_keepList)
head(locipassed_BREGMAJF_keepList)
locipassed_BREGMAJF_keepList <- drop.levels(locipassed_BREGMAJF_keepList)

################## Output list of Loci that passed MajAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/locipassed_BREGMAJF_keepList.txt", "wb")
write.table(locipassed_BREGMAJF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

######### make a joint list of the loci that passed both of the filters
length(locipassed_BREGMAF_keepList)
length(locipassed_BREGMAJF_keepList)
locipassed_MAF_MAJF_keepList <- union(locipassed_BREGMAF_keepList, locipassed_BREGMAJF_keepList )
length(locipassed_MAF_MAJF_keepList)


################## Output a list of the loci that passed either the MajAF or the MinAF filters############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/locipassed_MAF_MAJF_keepList.txt", "wb")
write.table(locipassed_MAF_MAJF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#########make a genotype set for the loci that passed either of the two filters so its easier to make their genepop files
chin_comb_RADGTseq_genotypes_MAF_MAJF_keep <- chin_comb_RADGTseq[,((colnames(chin_comb_RADGTseq) %in% locipassed_MAF_MAJF_keepList) | (colnames(chin_comb_RADGTseq) == "SampleName"))]
dim(chin_comb_RADGTseq_genotypes_MAF_MAJF_keep)
chin_comb_RADGTseq_genotypes_MAF_MAJF_keep[1:5,1:5]

################## Output the genotypes of the loci that passed either the MajAF or the MinAF filters############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADGTseq_genotypes_MAF_MAJF_keep.txt", "wb")
write.table(chin_comb_RADGTseq_genotypes_MAF_MAJF_keep,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)
