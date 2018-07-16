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

### Read in the data files we will need ########

#load the raw genepop file, edited for R. Includes all individuals and loci:  
chin_raw_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/batch_13.allCombined_NoCookInlet_R_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_raw_genepop)

#a population map that has all the individuals and the populations that they belong to. 
chin_all_popINFO <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/populationMap_NoCookInlet.txt", header =TRUE)
#colnames(chin_all_popINFO) <- c("SampleName", "Population", "VeryBroadRegion", "Location", "BroadRegion", "FineReportingGroup")
head(chin_all_popINFO)
dim(chin_all_popINFO)

#Garrett's HDPlot results for 19k loci
HDplot_GM_results <-read.delim("Z:/WORK/TARPEY/ChinookPanel/Chinook_HDplotResults.txt", header =TRUE)
head(HDplot_GM_results)
dim(HDplot_GM_results)

#import the snp positions for each locus
HD_plot_SNP_pos <- read.table("Z:/WORK/TARPEY/ChinookPanel/RAD_data/VCF_trueSNPpos.txt", header=FALSE, 
                              stringsAsFactors = FALSE, na.strings = "-" )
colnames(HD_plot_SNP_pos) <- c("POS", "Locus_ID", "Locus")
head(HD_plot_SNP_pos)

##Import the HDplot results for the remaining 3970 loci 
HDplot_new3970 <- read.table("Z:/WORK/TARPEY/ChinookPanel/newVCF/HDPlotResults_4.5.txt", header=TRUE, 
                              stringsAsFactors = FALSE, na.strings = "-" )
HDplot_new3970[1:5,1:10]


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
head(tags2exclude54)

#Replace the genotype with 0000 at any snp that is at a tag that has a SNP at position 54
dim(chin_raw_genepop)
dim(chin_filtering_genepop)

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
#rm(chin_raw_genepop)
#gc()

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
dim(lociMAF)

#if the loci does not pass the test, delete it
locipassedMAF<-vector()
locipassedMAF<-subset(lociMAF, lociMAF[,2]!=0)
dim(locipassedMAF)
head(locipassedMAF)

#what percent of the loci do we KEEP with this filter? 
Percent_locipassedMAF <- ((dim(locipassedMAF)[1])/(dim(lociMAF)[1]))*100
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
##retain loci that were at least 0.05 in any of the 4 broad regions
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
length(unique(locipassed_BREGMAF$Tag))

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
filtered_MAF_Genos_oneSNP_temp<-data.frame(keynames=names(filtered_MAF_Genos_oneSNP),value=filtered_MAF_Genos_oneSNP,row.names=NULL)
colnames(filtered_MAF_Genos_oneSNP_temp)<-c("Locus","MAF")
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
length(unique(oneMAF$Tag))
oneMAF<-as.data.frame(oneMAF)
oneMAF$Tag <- as.numeric(as.character(oneMAF$Tag))
oneMAF <- oneMAF[order(oneMAF$Tag, decreasing = FALSE),]
head(oneMAF)
dim(oneMAF)
#write.table(oneMAF,"Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/OneMAF.txt",quote=FALSE,row.names=FALSE)

#filter dataset to include only one SNP per tag loci
filtered_MAF_Genos_oneSNP<-filtered_MAF_Genos[,((colnames(filtered_MAF_Genos) %in% oneMAF$Locus) | (colnames(filtered_MAF_Genos) %in% colnames(chin_all_popINFO)))]
dim(filtered_MAF_Genos_oneSNP)

filtered_MAF_Genos[1:5,1:20]
dim(filtered_MAF_Genos)
filtered_MAF_Genos_oneSNP[1:5,1:11]

##################Get individual Genotype rate to filter individuals ###################

#get genotype rate per sample for filtered loci
filtered_MAF_Genos_oneSNP_inds<-apply(filtered_MAF_Genos_oneSNP,1,function(x) 1-(sum(x=="0000")/(dim(filtered_MAF_Genos_oneSNP)[2])))
filtered_MAF_Genos_oneSNP_inds<-data.frame(keyName=filtered_MAF_Genos_oneSNP$SampleName, value=filtered_MAF_Genos_oneSNP_inds, row.names=NULL)
colnames(filtered_MAF_Genos_oneSNP_inds)<-c("Sample","GenoRate")
dim(filtered_MAF_Genos_oneSNP_inds)
head(filtered_MAF_Genos_oneSNP_inds)

#plot ranked genotype rate for samples with filtered loci
filtered_MAF_Genos_oneSNP_inds_ranked<-filtered_MAF_Genos_oneSNP_inds[order(filtered_MAF_Genos_oneSNP_inds$GenoRate),]
filtered_MAF_Genos_oneSNP_inds_ranked$rank<-seq(1,dim(filtered_MAF_Genos_oneSNP_inds_ranked)[1],by=1)
ggplot()+geom_point(data=filtered_MAF_Genos_oneSNP_inds_ranked,aes(x=rank,y=GenoRate))+theme_bw() + 
  geom_hline(aes(yintercept=0.80),lty="dashed")+ggtitle("Sample Genotype Rate One SNP per Tag; Line Shows 80% Genotype Rate")
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/SampleGenotypeRateofLociGenotypedin80PCSamples.pdf")

#Keep samples with >=80% genotype rate
filteredSamples<-filtered_MAF_Genos_oneSNP_inds[filtered_MAF_Genos_oneSNP_inds$GenoRate>=0.80,]
filteredGenos_filteredSamples<-filtered_MAF_Genos_oneSNP[filtered_MAF_Genos_oneSNP$SampleName%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)
filteredGenos_filteredSamples[1:5,1:15]


############################Check the locus genotype rate again##################
#get genotype rate per locus for filtered sample dataset
locusGenoRate_filteredSamples<-apply(filteredGenos_filteredSamples[,c(7:dim(filteredGenos_filteredSamples)[2])],2,function(x) 1-(sum(x=="0000")/dim(filteredGenos_filteredSamples)[1]))
locusGenoRate_filteredSamples<-data.frame(keyName=names(locusGenoRate_filteredSamples), value=locusGenoRate_filteredSamples, row.names=NULL)
colnames(locusGenoRate_filteredSamples)<-c("Locus","GenoRate")
dim(locusGenoRate_filteredSamples)
head(locusGenoRate_filteredSamples)

#plot ranked genotype rate for samples with filtered loci
x_ranked<-locusGenoRate_filteredSamples[order(locusGenoRate_filteredSamples$GenoRate),]
x_ranked$rank<-seq(1,dim(x_ranked)[1],by=1)
ggplot()+geom_point(data=x_ranked,aes(x=rank,y=GenoRate))+ggtitle("Locus Genotype Rate of Filtered Samples")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/LocusGenotypeRateofFiltered465Samples.pdf")

## Should we re-filter the loci based on genotype rate now that there are indv missing? NO

###########Filter based on HDPLOT paralog status############  
##########these results are from Garrett Running HDPlot on his set of loci (19k) 
##So first I have to figure out what loci I am missing 

#merge the snp position with the HDplot results
head(HDplot_GM_results)
HDplot_GM_results_pos <- merge(HDplot_GM_results, HD_plot_SNP_pos, by= "Locus_ID") #this also sorts them by the locus ID
head(HDplot_GM_results_pos)

# ####WHich loci need to have HDplot run on a different SNP? 
# #
# #Compare my filtered Loci to Garretts HDplot results#
# loci_ct_filtered<-data.frame(str_split_fixed(locusGenoRate_filteredSamples$Locus,"_",2))
# names(loci_ct_filtered) <- c("Tag", "Snp")
# loci_ct_filtered$Locus <- locusGenoRate_filteredSamples$Locus
# head(loci_ct_filtered)
# 
# ###pull out loci that Garrett used
# HDplot_loci <- HDplot_GM_results_pos$Locus
# length(HDplot_loci)
# head(HDplot_loci)
# 
# #overlap between mine and his loci
# overlap_btw_HD_my <- intersect(HDplot_loci, loci_ct_filtered$Locus)
# length(overlap_btw_HD_my)
# head(overlap_btw_HD_my)
# 
# #In my set, not in his
# difference_btw_HD_my <- setdiff(loci_ct_filtered$Locus, HDplot_loci)
# length(difference_btw_HD_my)
# head(difference_btw_HD_my)
# 
# ####Write a table the loci that need to have a VCF file for them to run through HDPlot:
# outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/HDplot_SNPS_to_get_VCF.txt", "wb")
# write.table(difference_btw_HD_my,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)
# 
# #looking at the list of unique tags that I have in my 21055 set- to get all the SNPs from each of them to run through HDplot
# uniquetags <- locusGenoRate_filteredSamples$Locus
# uniquetags <- data.frame(str_split_fixed(locusGenoRate_filteredSamples$Locus,"_",2))
# uniquetags <- unique(uniquetags$X1)
# head(uniquetags)
# 
# ####Write a table of the unique tags that we have in the 21005 data set:
# outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/HDplot_tags_to_get_VCF.txt", "wb")
# write.table(uniquetags,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

### Combine the New HDplot Results of the missing 3970 loci with the original GM results 
#combine the two sets of : 
head(HDplot_new3970)
head(HDplot_GM_results_pos)
dim(HDplot_GM_results_pos)

#reorder the HDplot_GM_results_pos so that they line up with the ones for HDplot_new3970
HDplot_GM_results_pos_reorder <- HDplot_GM_results_pos[,c(12,2:10)]
colnames(HDplot_GM_results_pos_reorder) <- c("Locus_ID","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z","paralog")
head(HDplot_GM_results_pos_reorder)

HDPlot_combined_results <- rbind(HDplot_GM_results_pos_reorder, HDplot_new3970)
head(HDPlot_combined_results)
HDPlot_combined_results[19376:19382, 1:10]
dim(HDPlot_combined_results)

#make a list of the loci that are paralogs:
#looking at the column paralog- which are the yes paralogs? 
paralogs<-HDPlot_combined_results[HDPlot_combined_results$paralog != 0,]
dim(paralogs)

#filter the filtered data for paralogs

###paralog
HDplot_Paralogs <- paralogs$Locus_ID
length(HDplot_Paralogs)
head(HDplot_Paralogs)

#sort the list of paralogs numerically by tag
sorted_paralogs<-data.frame(str_split_fixed(HDplot_Paralogs,"_",2))
head(sorted_paralogs)
colnames(sorted_paralogs) <- c("tag", "pos")
sorted_paralogs$tag <- as.numeric(as.character(sorted_paralogs$tag))
sorted_paralogs <- sorted_paralogs[order(sorted_paralogs$tag),]
sorted_paralogs$locus <- paste(sorted_paralogs$tag, sorted_paralogs$pos, sep="_")
head(sorted_paralogs)
sorted_paralogs_locus <- sorted_paralogs$locus

#Keep samples that are not paralogs
filteredGenos_filteredSamples[1:5,1:15]
dim(filteredGenos_filteredSamples)
filteredGenos_filteredSamples_filtpara<-filteredGenos_filteredSamples[,(!(colnames(filteredGenos_filteredSamples) %in% HDplot_Paralogs))]

filteredGenos_filteredSamples<-filtered_MAF_Genos_oneSNP[filtered_MAF_Genos_oneSNP$SampleName%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples_filtpara)
filteredGenos_filteredSamples_filtpara[1:10,1:10]



########################################### Write the final genotypes and individuals to a file #########

####Write a table of the genotypes/individuals that have passed all the filters so far:
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/filteredGenos_filteredSamples_filtpara.txt", "wb")
write.table(filteredGenos_filteredSamples_filtpara,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Write a list of the individuals that have passed all the filters so far:
filteredSample_names <- filteredSamples$Sample
droplevels(filteredSample_names)
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/filteredSample_names.txt", "wb")
write.table(filteredSample_names,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

############# Write a txt file with the list of loci- use as whitelist or subset genepop file
#subset genepop file (for use with subset_genepop_by_SNPs.pl):
filtered_loci <- colnames(filteredGenos_filteredSamples_filtpara)
length(filtered_loci)
kept_filtered_loci <- filtered_loci[-c(1:6)]

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/kept_filtered_loci.txt", "wb")
write.table(filtered_loci,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

