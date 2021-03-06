### Running HDPlot on the chinook RAD data- all populations except Cook Inlet
###   This takes a VCF file and finds the true snp position, then finds the difference
###   between the set of loci we have the ones we are missing
###   A second VCF file for those loci is filtered to remove individuals in Cook Inlet
###   Then runs Garrett's code for HDplot to indentify paralogs
###  Garrett McKinney (and Carolyn Tarpey) | July 2018
### ---------------------------------------

#This includes R code written by Garrett McKinney to identify paralogs: HDPlot, the fast version
#It requires a VCF file, the results can be filtered for loci that we are interested later

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#input the raw VCF files, that have all the individuals
raw_vcf <- read.table("Z:/WORK/TARPEY/ChinookPanel/newVCF/batch_13_Tarpey_combined_R.vcf", header = TRUE)
head(raw_vcf)
dim(raw_vcf)

##input a list of the individuals and what populations they belong to 
all_inds <- read.table("Z:/WORK/TARPEY/ChinookPanel/newVCF/PopMap_allInds.txt", header = TRUE)
head(all_inds)

### filter out the individuals that we dont need from the VCF files to be run through HDplot
Inds_to_keep<-all_inds[all_inds$V.BroadRegion !="CookInlet",]
dim(Inds_to_keep)
intermediate_VCF <- raw_vcf[,colnames(raw_vcf)%in%Inds_to_keep$SampleName]
dim(intermediate_VCF)
intermediate_VCF[1:15,1:15]
refined_VCF <- cbind(raw_vcf[,c(1:9)],intermediate_VCF)
refined_VCF[1:15,1:15]

##export the refined VCF file, with the COOk Inlet samples removed
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/Refined_VCF_noCI.vcf", "wb")
write.table(refined_VCF,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE, sep="\t", eol="\n")
close(outputFile)

##export the refined VCF file, with the COOk Inlet samples removed to the DESKTOP
outputFile<-file("C:/Users/Carolyn/Desktop/Refined_VFC_noCI.vcf", "wb")
write.table(refined_VCF,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t", eol="\n")
close(outputFile)

#imput the combined VCF file, this will not load from the network drive- must be on the desktop
#you have to add the headers back to it too- they start with ##, and make sure it is tab delimited
vcfInput<-read.vcfR("C:/Users/Carolyn/Desktop/batch_13_Tarpey_combined.vcf")

vcfInput@fix[1:5,1:5]
vcfInput@gt[1:5,1:5]

HDPlot<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=9))
  colnames(HDplotTable)<-c("Locus_ID","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  #format allele reads from vcf data into matrix of comma separated values
  genos<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",3)[,1])
  rownames(genos)<-vcfData@fix[,3]
  reads<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",4)[,3])
  rownames(reads)<-vcfData@fix[,3] 
  #replace . with 0
  reads<-gsub("\\.","0",reads)
  alleleReads_1<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,1])
  alleleReads_2<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,2])
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  rownames(alleleReads_1)<-vcfData@fix[,3]
  rownames(alleleReads_2)<-vcfData@fix[,3]
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-genos
  hetMatrix<-apply(hetMatrix,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygotes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]
  
  #assign results to HDplotTable
  HDplotTable$Locus_ID<-vcfData@fix[,3]
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$het_perc<-hetPerc
  HDplotTable$std<-std
  HDplotTable$z<-z
  return(HDplotTable)
}

HDPlotResults<-HDPlot(vcfInput)
dim(HDPlotResults)

##plot the results of HDPlot in black and white
ggplot()+geom_point(data=HDPlotResults,aes(x=het_perc,y=z),alpha=0.15)+theme_bw()
ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/newVCF/BlackWhiteHDPlotHet_per_Z.pdf")

ggplot()+geom_point(data=HDPlotResults,aes(x=het_perc,y=ratio),alpha=0.15)+theme_bw()
ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/newVCF/BlackWhiteHDPlotHet_per_Ratio.pdf")

#set thresholds for marker types, these are from the Report that Garrett wrote to ADFG
thresh_H<-0.52
thresh_H_divDup<-0.85
thresh_D<-4.5
thresh_Dneg<--4.5

#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  H<-as.numeric(data[7])
  D<-as.numeric(data[9])
  if(is.na(D)){
    paralog<-NA
  }else if(H>=thresh_H_divDup){
    paralog<-2
  }else if(H>=thresh_H|D>thresh_D|D<thresh_Dneg){
    paralog<-1
  }else{
    paralog<-0
  }
  return(paralog)
}

HDPlotResults$paralog<-apply(HDPlotResults,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup)
HDPlotResults$paralog<-as.factor(HDPlotResults$paralog)


ggplot() + geom_point(data=HDPlotResults, aes(x=het_perc,y=z,color=paralog), alpha=0.15) + scale_color_manual(values=c("blue","red","dark green"), 
  labels = c("Singletons", "Duplicates", "Diverged Duplicates","")) + theme_bw() + theme(legend.title=element_blank(),
  legend.text = element_text(size= 15))+ guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/newVCF/ColorHDPlotHet_per_Z_4.5.pdf") 

ggplot()+geom_point(data=HDPlotResults,aes(x=het_perc,y=ratio,color=paralog),alpha=0.15)+scale_color_manual(values=c("blue","red","dark green"),
  labels = c("Singletons", "Duplicates", "Diverged Duplicates","")) + theme_bw() + theme(legend.title=element_blank(),
  legend.text = element_text(size= 15))+ guides(colour = guide_legend(override.aes = list(alpha = 1)))
ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/newVCF/ColorHDPlotHet_per_Ratio4.5.pdf")  

#count the number of each type of locus                                                                                                                                                                                                                                 
sum(HDPlotResults$paralog==0,na.rm=TRUE) # singletons
sum(HDPlotResults$paralog==1,na.rm=TRUE) # duplicates
sum(HDPlotResults$paralog==2,na.rm=TRUE) # diverged duplicates

# ##convert the position number in the vcf file to the accurate position of the snp in the tag.
# ### ours are already included
# head(HDPlotResults)
# vcfInput@fix[1:5,1:5]
# wrongPositions<-as.data.frame(vcfInput@fix)
# wrongPositions$POS<-as.numeric(as.character(wrongPositions$POS))
# correctPositions<-(wrongPositions$POS-2)%%94
# correctPositions
# 
# #put these positions in a data frame
# VCF_paralogs<-data.frame(value=correctPositions,row.names=names(correctPositions))
# colnames(VCF_paralogs)<-c("Position")
# 
# #concatenate the locus names to the results.
# wrongPositions$ID<-as.numeric(as.character(wrongPositions$ID))
# Tag <- wrongPositions$ID
# length(Tag)
# VCF_paralogs$Tag<-Tag
# head(VCF_paralogs)
# 
# #add a column that is the full locus name
# VCF_paralogs$Locus <- paste(VCF_paralogs$Tag,VCF_paralogs$Position,sep="_")
# 
# #Add the paralog identity to the data frame
# VCF_paralogs$Identity<-HDPlotResults$paralog
# head(VCF_paralogs)


VCF_paralogs <- HDPlotResults
head(VCF_paralogs)

#export a table of the HDPLOT results
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/HDPlotResults_4.5.txt", "wb")
write.table(HDPlotResults,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)


####THe following is optional ####

#add an additional column with the identity in words
VCF_paralogs$Identity_2<-as.character(HDPlotResults$paralog)
VCF_paralogs$Identity_2[VCF_paralogs$Identity_2==0] <- "singleton"
VCF_paralogs$Identity_2[VCF_paralogs$Identity_2==1] <- "duplicate"
VCF_paralogs$Identity_2[VCF_paralogs$Identity_2==2] <- "div_duplicate"
head(VCF_paralogs)

#export a table of the loci and their paralog status
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/VCF_paralogs_4.5.txt", "wb")
write.table(VCF_paralogs,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

#sort the table and create lists of each of the loci for each paralog status
singletons <- VCF_paralogs[VCF_paralogs$paralog==0,]
singletons <- as.data.frame(singletons)
singletons <- na.omit(singletons)

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/singletons_4.5.txt", "wb")
write.table(singletons,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

duplicates <- VCF_paralogs[VCF_paralogs$paralog==1,]
duplicates <- as.data.frame(duplicates)
duplicates <- na.omit(duplicates)

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/duplicates_4.5.txt", "wb")
write.table(duplicates,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

div_duplicates <- VCF_paralogs[VCF_paralogs$paralog==2,]
div_duplicates <- as.data.frame(div_duplicates)
div_duplicates <- na.omit(div_duplicates)

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/div_duplicates_4.5.txt", "wb")
write.table(div_duplicates,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)


##Any paralog
any_paralog <- VCF_paralogs[VCF_paralogs$paralog!=0,]
any_paralog <- as.data.frame(any_paralog)
any_paralog <- na.omit(any_paralog)
dim(any_paralog)
head(any_paralog)

outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/newVCF/any_paralog_4.5.txt", "wb")
write.table(any_paralog,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)



