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
chin_raw_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/batch_13.allCombined_NoCookInlet_R_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_raw_genepop)

#a population map that has all the individuals and the populations that they belong to. 
chin_all_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/populationMap_NoCookInlet.txt", header =TRUE)
#colnames(chin_all_popINFO) <- c("SampleName", "Population", "VeryBroadRegion", "Location", "BroadRegion", "FineReportingGroup")
head(chin_all_popINFO)
dim(chin_all_popINFO)

######################################
#remove the commas from the sample names, and the X from the column names
chin_raw_genepop[1:15,1:15]
rownames(chin_raw_genepop)<-gsub(",","",rownames(chin_raw_genepop))
colnames(chin_raw_genepop)<-gsub("X","",colnames(chin_raw_genepop))
chin_raw_genepop[1:15,1:15]

####Remove Tags with SNPS at position 54#####
#due to sequencing error

#get unique tags from full dataset
allLoci<-data.frame(str_split_fixed(colnames(chin_raw_genepop),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(chin_raw_genepop)
length(allLoci$Locus)
length(unique(allLoci$Tag))
dim(allLoci)
allLoci[1:15,1:3]

#which tags have a SNP at Position 54?
tags_pos_54 <- allLoci[allLoci$SNP == 54,]
dim(tags_pos_54)
tags2exclude54 <- unique(tags_pos_54$Tag)
length(tags2exclude54)
head(tags_pos_54)

#Replace the genotype with 0000 at any snp that is at a tag that has a SNP at position 54
dim(chin_raw_genepop)

chin_filtering_genepop <-chin_raw_genepop

for (a in 1:length(allLoci$Locus)) {
  temp1 <- (str_split_fixed(colnames(chin_filtering_genepop)[a],"_",2))
  if (temp1[1] %in% tags2exclude54){
    chin_filtering_genepop[,a]<-"0000"
  } else {
    next()
  }
}

#the above bogged down my machine so I deleted the raw genotypes and cleared the memory
# rm(chin_raw_genepop)
# gc()

####### genotype rate per locus################

locusGenoRate<-apply(chin_filtering_genepop,2,function(x) 1-(sum(x=="0000")/dim(chin_filtering_genepop)[1]))
locusGenoRate<-data.frame(keyName=names(locusGenoRate), value=locusGenoRate, row.names=NULL)
colnames(locusGenoRate)<-c("Locus","GenoRate")
dim(locusGenoRate)
head(locusGenoRate)
#write.table(locusGenoRate,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/LocusGenoRate.txt",quote=FALSE,row.names=FALSE)

#filterLoci with >= 80% genotype rate
filteredLoci<-locusGenoRate[locusGenoRate$GenoRate>=0.80,]
dim(filteredLoci)
head(filteredLoci)

#get unique tags from 80% genoRate filtered dataset
retainedLoci_80perc<-data.frame(str_split_fixed(filteredLoci$Locus,"_",2))
colnames(retainedLoci_80perc)<-c("Tag","SNP")
retainedLoci_80perc$Locus<-filteredLoci$Locus
head(retainedLoci_80perc)
length(unique(retainedLoci_80perc$Tag))

# #format dataset to remove X from locus names
# filteredLociIDs<-filteredLoci$Locus
# filteredLociIDs<-gsub("X","",filteredLociIDs)

#filter dataset to include only filtered loci
filteredGenos<-chin_filtering_genepop[,colnames(chin_filtering_genepop)%in%filteredLoci$Locus]
dim(filteredGenos)

chin_filtering_genepop[1:5,1:5]
filteredGenos[1:5,1:5]
dim(filteredGenos)

#write.table(filteredGenos,"Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/FirstFiltering/filteredGenos_just80PCgenorate_.txt",quote=FALSE,row.names=TRUE)

##############Minor Allele Frequency MAF#############

#assign each sample to a population 
SampleName <- rownames(filteredGenos)
filteredGenos <- cbind(SampleName, filteredGenos)
#filteredGenos[,2:7] <- NULL
filteredGenos<-merge(chin_all_popINFO,filteredGenos, by="SampleName")
filteredGenos[1:5,1:15]
dim(filteredGenos)


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

######### Test MAF >= 0.05 in at least one of 16 pops ############

#create a dataframe for our MAF results
npops<- length(unique(filteredGenos$Population))
nloci<- dim(retainedLoci_80perc)[1]

filteredGenos_POP_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(filteredGenos_POP_MAF)<-as.vector(unique(filteredGenos$Population)) #name the rows by the population names
colnames(filteredGenos_POP_MAF)<-retainedLoci_80perc$Locus
filteredGenos_POP_MAF[1:6,1:6]
dim(filteredGenos_POP_MAF)


##### Get the MAF for each locus by population 
for (i in 1:16){ #16 is the unique number of populations we have
  tempPop<-(unique(filteredGenos$Population)[i])
  tempGeno<-subset(filteredGenos, Population == tempPop)
  # tempGeno[,1:15]
  # dim(tempGeno)
  tempMAF<-apply(tempGeno[,c(7:dim(filteredGenos)[2])],2,calculateMAF)
  filteredGenos_POP_MAF[i,]<-c(tempMAF)
}

filteredGenos_POP_MAF[1:16,1:5]
dim(filteredGenos_POP_MAF)
#colnames(filteredGenos_POP_MAF[,1])

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 16 populations
lociMAF<-vector()
lociMAF<-apply(filteredGenos_POP_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
lociMAF<-data.frame(keynames=names(lociMAF), value=lociMAF, row.names = NULL)
colnames(lociMAF)<-c("Locus", "PopsMAF")

#if the loci does not pass the test, delete it
locipassedMAF<-vector()
locipassedMAF<-subset(lociMAF, lociMAF[,2]!=0)
dim(locipassedMAF)
head(locipassedMAF)

#what percent of the loci do we KEEP with this filter? 
Percent_locipassedMAF <- ((dim(locipassedMAF)[1])/70530)*100
Percent_locipassedMAF

#what percent of the tags do we KEEP with this filter? 
locipassedMAF<-data.frame(str_split_fixed(locipassedMAF$Locus,"_",2))
colnames(locipassedMAF)<-c("Tag","SNP")
UNIQUEtags_locipassedMAF<-unique(locipassedMAF$Tag)
length(UNIQUEtags_locipassedMAF)
Percent_UNIQUEtags_locipassedMAF <- ((length(UNIQUEtags_locipassedMAF))/40924)*100
Percent_UNIQUEtags_locipassedMAF

######### Test MAF >= 0.05 in at least one of the 4 main regions: BroadReportingGroup (BroadRegion, column 5) ############
colnames(filteredGenos[,c(1:10)])

#create a dataframe for our MAF results
npops<- length(unique(filteredGenos$BroadRegion))
nloci<- dim(retainedLoci_80perc)[1]

filteredGenos_BREG_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(filteredGenos_BREG_MAF)<-as.vector(unique(filteredGenos$BroadRegion)) #name the rows by the population names
colnames(filteredGenos_BREG_MAF)<-retainedLoci_80perc$Locus
dim(filteredGenos_BREG_MAF)
filteredGenos_BREG_MAF[1:4,1:6]

##### Get the MAF for each locus by region 
for (i in 1:4){ #4 is the unique number of Broad Reporting Groups we have
  tempRegion<-(unique(filteredGenos$BroadRegion)[i])
  tempGeno<-subset(filteredGenos, BroadRegion == tempRegion)
  tempMAF<-apply(tempGeno[,c(7:dim(filteredGenos)[2])],2,calculateMAF)
  filteredGenos_BREG_MAF[i,]<-c(tempMAF)
}

filteredGenos_BREG_MAF[1:4,1:5]
dim(filteredGenos_BREG_MAF)
#colnames(filteredGenos_BREG_MAF[,1])

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 16 populations
loci_BREGMAF<-vector()
loci_BREGMAF<-apply(filteredGenos_BREG_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
loci_BREGMAF<-data.frame(keynames=names(loci_BREGMAF), value=loci_BREGMAF, row.names = NULL)
colnames(loci_BREGMAF)<-c("Locus", "BroadRegionzMAF")
head(loci_BREGMAF)

#if the loci does not pass the test, delete it
locipassed_BREGMAF<-vector()
locipassed_BREGMAF<-subset(loci_BREGMAF, loci_BREGMAF[,2]!=0)
dim(locipassed_BREGMAF)
head(locipassed_BREGMAF)

## what percent of the loci we keep with this filter? 
Percent_locipassed_BREGMAF <- ((dim(locipassed_BREGMAF)[1])/70530)*100
Percent_locipassed_BREGMAF

#get list of loci that made it past this filter
locipassed_BREGMAF_keepList <- locipassed_BREGMAF$Locus
length(locipassed_BREGMAF_keepList)
head(locipassed_BREGMAF_keepList)

#get unique tags from initial filtered dataset
locipassed_BREGMAF<-data.frame(str_split_fixed(locipassed_BREGMAF$Locus,"_",2))
colnames(locipassed_BREGMAF)<-c("Tag","SNP")

## what percent of the loci we keep with this filter? 
UNIQUEtags_locipassed_BREGMAF<-unique(locipassed_BREGMAF$Tag)
length(UNIQUEtags_locipassed_BREGMAF)
Percent_UNIQUEtags_locipassed_BREGMAF <- ((length(UNIQUEtags_locipassed_BREGMAF))/40924)*100
Percent_UNIQUEtags_locipassed_BREGMAF

################# Filter out Loci that do not meet the MAF threshold of 0.05 in 1 or more regions ##########################
filtered_MAF_Genos<-filteredGenos[,((colnames(filteredGenos) %in% locipassed_BREGMAF_keepList) | (colnames(filteredGenos) %in% colnames(chin_all_popINFO)))]
dim(filtered_MAF_Genos)
filtered_MAF_Genos[1:4,1:15]

################## Output lists of Loci that passed MAF filter############
#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/locipassed_BroadREGMAF_keepList.txt", "wb")
write.table(locipassed_BREGMAF_keepList,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


################## Create a dataset that is one SNP per tag- Highest MAF ############
# use the genotypes from filtered_MAF_Genos

#the function to calculate MAF
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

filtered_MAF_Genos[1:5,1:15]

#this calculates the MAF overall. 
filtered_MAF_Genos_oneSNP<-apply(filtered_MAF_Genos[,c(7:(dim(filtered_MAF_Genos)[2]))],2,calculateMAF)
head(filtered_MAF_Genos_oneSNP)

#put the results in a dataframe
filtered_MAF_Genos_oneSNP_temp<-data.frame(value=filtered_MAF_Genos_oneSNP,row.names=names(filtered_MAF_Genos_oneSNP))
colnames(filtered_MAF_Genos_oneSNP_temp)<-c("MAF")
head(filtered_MAF_Genos_oneSNP_temp)

#concatenate the locus names to the results. 
Loci_temp<-colnames(filtered_MAF_Genos[,c(7:(dim(filtered_MAF_Genos)[2]))])
length(Loci_temp)
filtered_MAF_Genos_oneSNP_temp$Locus<-Loci_temp
head(filtered_MAF_Genos_oneSNP_temp)

#split the tags and snp positions
filtered_MAF_Genos_oneSNP_temptags<-data.frame(str_split_fixed(filtered_MAF_Genos_oneSNP_temp$Locus,"_",2))
colnames(filtered_MAF_Genos_oneSNP_temptags)<-c("Tag","SNP")
head(filtered_MAF_Genos_oneSNP_temptags)
filtered_MAF_Genos_oneSNP_temp$Tag<-filtered_MAF_Genos_oneSNP_temptags$Tag
filtered_MAF_Genos_oneSNP_temp$SNP<-filtered_MAF_Genos_oneSNP_temptags$SNP
head(filtered_MAF_Genos_oneSNP_temp)
dim(filtered_MAF_Genos_oneSNP_temp)

#write.table(filtered_MAF_Genos_oneSNP_temp,"Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/filteredforBroadRegMAF_MAFglobal_Genos.txt",quote=FALSE,row.names=FALSE)

#Retain the SNP with the highest MAF per tag
# This one gives an output with the correct number of unique tags, and spot checked in excel
oneMAF<- filtered_MAF_Genos_oneSNP_temp %>% group_by(Tag) %>% slice(which.max(MAF))
head(oneMAF)
oneMAF<-as.data.frame(oneMAF)
head(oneMAF)
dim(oneMAF)

#write.table(oneMAF,"Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/OneMAF.txt",quote=FALSE,row.names=FALSE)

#filter dataset to include only one SNP per tag loci
length(oneMAF$Locus)

filtered_MAF_Genos_oneSNP<-filtered_MAF_Genos[,oneMAF$Locus]
dim(filtered_MAF_Genos_oneSNP)

filtered_MAF_Genos[1:5,7:15]
filtered_MAF_Genos_oneSNP[1:5,7:15]

####<------------------------------------------------START HERE,
############# Write a txt file with the list of loci- use as whitelist or subset genepop file

#subset genepop file (for use with subset_genepop_by_SNPs.pl):

#format dataset to remove X from locus names
oneMAF_loci<-gsub("X","",oneMAF$Locus)
head(oneMAF_loci)

outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/oneSNPpertag_subsetgenepop.txt", "wb")
write.table(oneMAF_loci,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


##################### Run Genepop with the genepop file for HWE and Fis and import the results here for filtering based on those 

#### Take the original Genepop file that was formatted for import at the beginning of this code and make it on snp per line
###  Then subset it with the list that was exported above, and run it in Genepop for HWE
###  Take the output of that and use the perl script get_HWresults_from_genepop.pl to get the results from the INF output file

HWE_table_Pre_Filtered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/batch_4_31485LOCI_HWE_results.txt" ,
                                     stringsAsFactors = FALSE, row.names= 1)
HWE_pops <-c("AMUR_10", "AMUR_11", "SUSIT_13", "HAYLY_09", "HAYLY_10", "KOPE_91", "KOPE_96", "KUSHI_06", "KUSHI_07",
             "LAKEL_06", "LAKEL_07", "NOME_91", "NOME_94", "SNOH_03", "SNOH_96", "SUSIT_14", "TAUY_09", "TAUY_12")

colnames(HWE_table_Pre_Filtered)<-c("AMUR_10", "AMUR_11", "SUSIT_13", "HAYLY_09", "HAYLY_10", "KOPE_91", "KOPE_96", "KUSHI_06", "KUSHI_07",
                                    "LAKEL_06", "LAKEL_07", "NOME_91", "NOME_94", "SNOH_03", "SNOH_96", "SUSIT_14", "TAUY_09", "TAUY_12")
head(HWE_table_Pre_Filtered)
dim(HWE_table_Pre_Filtered)
nrow(HWE_table_Pre_Filtered)

##### Filter the loci by the HWE
##retain loci that were at least 0.05 in at least 9 of the populations

loci_HWE_blank_test<-vector()
loci_HWE_blank_test<-apply(HWE_table,1,function(x) sum(x <=0.05, na.rm=TRUE)-sum(x =="-", na.rm=TRUE))
loci_HWE_blank_test<-data.frame(keynames=names(loci_HWE_blank_test), value=loci_HWE_blank_test, row.names = NULL)
colnames(loci_HWE_blank_test)<-c("Locus", "failedHWE-blanks")
head(loci_HWE_blank_test)
dim(loci_HWE_blank_test)

#if the loci does not pass the test, delete it from this group
locipassedHWE_test<-vector()
locipassedHWE_test<-subset(loci_HWE_blank_test, loci_HWE_blank_test[,2]<=9)
dim(locipassedHWE_test)
head(locipassedHWE_test)

# #if the loci does not pass the test put it into a list for later
# lociFailedHWE_test<-vector()
# lociFailedHWE_test<-subset(loci_HWE_blank_test, loci_HWE_blank_test[,2]>=9)
# dim(lociFailedHWE_test)
# head(lociFailedHWE_test)
# outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/Loci_FAILED_HWE.txt", "wb")
# write.table(lociFailedHWE_test[,1],outputFile,quote=FALSE,row.names=TRUE,col.names=FALSE,eol="\n")
# close(outputFile)

# #Retain the Loci that Passed HWE
# HWE_table_Post_Filtered <-  HWE_table_Pre_Filtered[row.names(HWE_table_Pre_Filtered)%in%locipassedHWE_test$Locus, ] #this subsets the matrix by the row names in the list
# dim(HWE_table_Post_Filtered)
# head(HWE_table_Post_Filtered)


#####################Test to see that what we got in One SNP per tag by MAF has the 16681 loci that we want

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci names in oneMAF_temp
locipassedHWE_test_comp<- locipassedHWE_test$Locus
head(locipassedHWE_test_comp)
length(locipassedHWE_test_comp)

#loci that are in 14637 and 30088, the loci that passed HWE
length(in_16662_and_31485)
in_14637_and_30088 <- intersect(in_16662_and_31485, locipassedHWE_test_comp)
length(in_14637_and_30088)

diff_14637_and_30088 <- setdiff(in_16662_and_31485, locipassedHWE_test_comp)
length(diff_14637_and_30088)

###########################################

##################subsample the genotypes based on the Loci that passed HWE
#all_newGenos <- temp_allnewGenos #reset

#format genotype dataset to remove X from locus names
#temp_allnewGenos <- all_newGenos
colnames(all_newGenos)<-gsub("X","",colnames(all_newGenos))
all_newGenos[1:5,1:5]

#filter dataset to include only loci that passed HWE
all_newGenos_HWE_filtered<-all_newGenos[,colnames(all_newGenos)%in%locipassedHWE_test$Locus]
dim(all_newGenos_HWE_filtered)
all_newGenos_HWE_filtered[1:5,1:5]


####Write list of loci that passed all filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/Loci_passed_HWE_filtered.txt", "wb")
write.table(colnames(all_newGenos_HWE_filtered),outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



#get genotype rate per sample for filtered loci
all_newGenos_HWE_inds<-apply(all_newGenos_HWE_filtered,1,function(x) 1-(sum(x=="0000")/dim(all_newGenos_HWE_filtered)[2]))
all_newGenos_HWE_inds<-data.frame(keyName=names(all_newGenos_HWE_inds), value=all_newGenos_HWE_inds, row.names=NULL)
colnames(all_newGenos_HWE_inds)<-c("Sample","GenoRate")
dim(all_newGenos_HWE_inds)
head(all_newGenos_HWE_inds)


#plot ranked genotype rate for samples with filtered loci
all_newGenos_HWE_inds_ranked<-all_newGenos_HWE_inds[order(all_newGenos_HWE_inds$GenoRate),]
all_newGenos_HWE_inds_ranked$rank<-seq(1,dim(all_newGenos_HWE_inds_ranked)[1],by=1)
ggplot()+geom_point(data=all_newGenos_HWE_inds_ranked,aes(x=rank,y=GenoRate))+theme_bw() + 
  geom_hline(aes(yintercept=0.80),lty="dashed")+ggtitle("Sample Genotype Rate One SNP per Tag; Line Shows 80% Genotype Rate")
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/SampleGenotypeRateofLociGenotypedin80PCSamples.pdf")


#Keep samples with >=80% genotype rate
filteredSamples<-all_newGenos_HWE_inds[all_newGenos_HWE_inds$GenoRate>=0.80,]
filteredGenos_filteredSamples<-all_newGenos_HWE_filtered[rownames(all_newGenos_HWE_filtered)%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)


#get genotype rate per locus for filtered sample dataset
locusGenoRate_filteredSamples<-apply(filteredGenos_filteredSamples,2,function(x) 1-(sum(x=="0000")/dim(all_newGenos_HWE_inds)[1]))
locusGenoRate_filteredSamples<-data.frame(keyName=names(locusGenoRate_filteredSamples), value=locusGenoRate_filteredSamples, row.names=NULL)
colnames(locusGenoRate_filteredSamples)<-c("Locus","GenoRate")
dim(locusGenoRate_filteredSamples)
head(locusGenoRate_filteredSamples)

####Write a table of the genotypes/individuals that have passed all the filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/FilteredGenos_FilteredInds.txt", "wb")
write.table(filteredGenos_filteredSamples,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Write a list of the individuals that have passed all the filters so far:
filteredSample_names <- row.names(filteredGenos_filteredSamples)
filteredSample_names <- as.list(filteredSamples$Sample)
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/filteredSample_names.txt", "wb")
write.table(filteredSample_names,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

####Write a list of the individuals that have FAILED the filters so far:
#samples with >=80% genotype rate
FAILED_filteredSamples<-all_newGenos_HWE_inds[all_newGenos_HWE_inds$GenoRate<0.80,]
dim(FAILED_filteredSamples)
head(FAILED_filteredSamples)
FAILED_filteredSample_names <-as.list(FAILED_filteredSamples$Sample)
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/FAILED_filteredSamples.txt", "wb")
write.table(FAILED_filteredSample_names,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


######## Should we re-filter the loci based on genotype rate now that there are indv missing? NO

#plot ranked genotype rate for samples with filtered loci
x_ranked<-locusGenoRate_filteredSamples[order(locusGenoRate_filteredSamples$GenoRate),]
x_ranked$rank<-seq(1,dim(x_ranked)[1],by=1)
ggplot()+geom_point(data=x_ranked,aes(x=rank,y=GenoRate))+ggtitle("Locus Genotype Rate of Filtered 465 Samples")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/LocusGenotypeRateofFiltered465Samples.pdf")

#Which loci are genotyped at less than 80% Genotype rate? 
x_80<-locusGenoRate_filteredSamples[locusGenoRate_filteredSamples$GenoRate<0.80,]
x_80_sort <-x_80[order(x_80$GenoRate, decreasing = FALSE),]
dim(x_80_sort)
head(x_80_sort)






