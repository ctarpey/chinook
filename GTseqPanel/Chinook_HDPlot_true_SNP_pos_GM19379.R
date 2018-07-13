### Finding the true SNP pos from a VCF file: chinook GTseq raw rad loci
###   This takes a VCF file  
###   
###  Garrett McKinney (and Carolyn Tarpey) | July 2018
### ---------------------------------------

#This includes R code written by Garrett McKinney 
#It requires a VCF file, the results can be filtered for loci that we are interested in later

library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#imput the combined VCF file, this will not load from the network drive- must be on the desktop
vcfInput<-read.vcfR("C:/Users/Carolyn/Desktop/batch_13_full_WAKreport.vcf")

vcfInput@fix[1:5,1:5]

##convert the position number in the vcf file to the accurate position of the snp in the tag. 
head(HDPlotResults)
vcfInput@fix[1:5,1:5]
wrongPositions<-as.data.frame(vcfInput@fix)
wrongPositions$POS<-as.numeric(as.character(wrongPositions$POS))
head(wrongPositions)
correctPositions<-(wrongPositions$POS-2)%%94 
correctPositions

#put these positions in a data frame 
VCF_paralogs<-data.frame(value=correctPositions,row.names=names(correctPositions))
colnames(VCF_paralogs)<-c("Position")
head(VCF_paralogs)

#concatenate the locus names to the results. 
wrongPositions$ID<-as.numeric(as.character(wrongPositions$ID))
Tag <- wrongPositions$ID
length(Tag)
VCF_paralogs$Tag<-Tag
head(VCF_paralogs)

#add a column that is the full locus name 
VCF_paralogs$Locus <- paste(VCF_paralogs$Tag,VCF_paralogs$Position,sep="_")

#export a table of the loci and their paralog status
outputFile<-file("Z:/WORK/TARPEY/ChinookPanel/RAD_data/VCF_trueSNPpos.txt", "wb")
write.table(VCF_paralogs,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



