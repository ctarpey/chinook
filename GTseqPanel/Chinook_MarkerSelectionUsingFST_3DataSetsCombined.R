### To select markers for a refined Chinook panel 
###   Using Ranked FST of one snp per tag loci, combining the 3 data sets <--- Method 2
###    
### to identify markers that distinguish the Kusk from Nush
###
###   Carolyn Tarpey | August 2018
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
#They already have primers made, so they dont have to have the FST per locus-theyre good wrapped up like this

###### Overlapping loci
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
dim(allCombined_FST_file)
allCombined_FST_file[1:5,1:5]


## RADGTseq data set estimated FST
RADGTseq_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADGTSeq_nonduplicated_FST_results.txt", header=TRUE,  row.names=NULL,
                              stringsAsFactors = FALSE, na.strings = "-" )
dim(RADGTseq_FST_file)
RADGTseq_FST_file[1:5,1:5]


## GTseqONLY data set estimated FST
GTseqONLY_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/GTSeqONLY_nonduplicated_FST_results.txt", header=TRUE, row.names=NULL,
                               stringsAsFactors = FALSE, na.strings = "-" )
dim(GTseqONLY_FST_file)
GTseqONLY_FST_file[1:5,1:5]


#there are different populations included in the newRAD data set than in the RADGTSeq and the GTseqONLY sets. This doesnt matter for the 
#Kusk or the Nush individual or joint analyses, just for the All. I removed the Yukon and Norton Sound pops so that there are the same populations
#represented in the newRAD as there are in the RADGTseq, to make a cleaner comparison between the two. 
AllmatchRADGTseq_FST_file<-read.table("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/RADallCombined_nondup_FST_results_ALLmatchRADGTSeq.txt", header=TRUE, row.names=NULL,
                               stringsAsFactors = FALSE, na.strings = "-" )

dim(AllmatchRADGTseq_FST_file)
AllmatchRADGTseq_FST_file[1:5,1:2]

###################################### Marker selection based on ranked FSt of loci 
##look at the FST of the snps and choose the ones that have the highest FST for the primer design and panel development. 

######################################Identify which loci are in the allCombined that are also in the RADGTseq
##This will require converting the names of one of the sets of loci. 

RADGTseq_loci <- RADGTseq_FST_file$row.names
head(RADGTseq_loci)

# ### Write the RADGTseq loci to a table
# outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/RADGTseq_loci", "wb")
# write.table(RADGTseq_loci,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

#after I exported this list I took it to excel and separated out the loci with these naming schemes: RAD392 Ots_RAD22318 Ots_RAD6618-57 from the rest of the loci
#I got rid of the RAD, Ots_RAD and the -# to just leave the tag number. I made a text file with two columns, the original name and the tag number and imported that 
##here: 

RADGTseq_locus_tags <- read.table("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/RADGTseq_locus_tags_to_match_newRAD.txt" , header=TRUE)
dim(RADGTseq_locus_tags)
head(RADGTseq_locus_tags)

#separate out the tag names from the allCombined loci 
allCombined_loci <- allCombined_FST_file$row.names
allCombined_loci_tags<-data.frame(str_split_fixed(allCombined_loci,"_",2))
allCombined_loci_tags$X3<-allCombined_loci
colnames(allCombined_loci_tags) <- c("TAG", "SNP","Locus")
head(allCombined_loci_tags)
dim(allCombined_loci_tags)

RADGTseq_RADnew_locus_overlap <- intersect(allCombined_loci_tags$TAG,RADGTseq_locus_tags$TAG )
length(RADGTseq_RADnew_locus_overlap)
head(RADGTseq_RADnew_locus_overlap)

## remove the overlapping loci from the RADNEW datasets
## not the older set, because it doesnt capture the haplotypes the way the old one does. 

###Make a data frame of the ones to keep 

allCombined_loci_tags_to_keep <- allCombined_loci_tags[!(allCombined_loci_tags$TAG %in% RADGTseq_RADnew_locus_overlap),]
dim(allCombined_loci_tags_to_keep)
head(allCombined_loci_tags_to_keep)

allCombined_FST_tags_to_keep <- allCombined_FST_file[ allCombined_FST_file$row.names %in% allCombined_loci_tags_to_keep$Locus,]
dim(allCombined_FST_tags_to_keep)
head(allCombined_FST_tags_to_keep)

# ### Write the allCombined_loci_tags_to_keep to a file 
# outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/allCombined_loci_tags_to_keep.txt", "wb")
# write.table(allCombined_loci_tags_to_keep,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

######################################Identify which loci are in the allCombined that are also in the RADGTseq
##This will require converting the names of one of the sets of loci. 

AllmatchRADGTseq_loci <- AllmatchRADGTseq_FST_file$row.names
head(AllmatchRADGTseq_loci)

AllmatchRADGTseq_loci_tags<-data.frame(str_split_fixed(AllmatchRADGTseq_loci,"_",2))
AllmatchRADGTseq_loci_tags$X3<-AllmatchRADGTseq_loci
colnames(AllmatchRADGTseq_loci_tags) <- c("TAG", "SNP","Locus")
head(AllmatchRADGTseq_loci_tags)
dim(AllmatchRADGTseq_loci_tags)

AllmatchRADGTseq_RADnew_locus_overlap <- intersect(RADGTseq_locus_tags$TAG,AllmatchRADGTseq_loci_tags$TAG)
length(AllmatchRADGTseq_RADnew_locus_overlap)
head(AllmatchRADGTseq_RADnew_locus_overlap)

## remove the overlapping loci from the AllmatchRADGTseq_RADnew datasets
## not the older set, because the new doesnt capture the haplotypes the way the old one does. 

###Make a data frame of the ones to keep 

AllmatchRADGTseq_loci_tags_to_keep <- AllmatchRADGTseq_loci_tags[!(AllmatchRADGTseq_loci_tags$TAG %in% AllmatchRADGTseq_RADnew_locus_overlap),]
dim(AllmatchRADGTseq_loci_tags_to_keep)
head(AllmatchRADGTseq_loci_tags_to_keep)

AllmatchRADGTseq_FST_tags_to_keep <- AllmatchRADGTseq_FST_file[ AllmatchRADGTseq_FST_file$row.names %in% AllmatchRADGTseq_loci_tags_to_keep$Locus,]
dim(AllmatchRADGTseq_FST_tags_to_keep)
head(AllmatchRADGTseq_FST_tags_to_keep)


######################################Combine them into one big data set by group
##put the three different data sets together, by population group 

head(allCombined_FST_tags_to_keep)
dim(allCombined_FST_tags_to_keep)

head(GTseqONLY_FST_file)
dim(GTseqONLY_FST_file)

head(RADGTseq_FST_file)
dim(RADGTseq_FST_file)

#add a column to each that has their identity, so that later we know what set of loci they belong to.
allCombined_FST_tags_to_keep$Identity <- "NewRAD"
head(allCombined_FST_tags_to_keep)
dim(allCombined_FST_tags_to_keep)

GTseqONLY_FST_file$Identity <- "GTseqONLY"
head(GTseqONLY_FST_file)
dim(GTseqONLY_FST_file)

RADGTseq_FST_file$Identity <- "RADGTseq"
head(RADGTseq_FST_file)
dim(RADGTseq_FST_file)

intermediate_1 <- rbind(allCombined_FST_tags_to_keep, GTseqONLY_FST_file )
head(intermediate_1)
dim(intermediate_1)

FST_3Combined <- rbind(intermediate_1, RADGTseq_FST_file )
head(FST_3Combined)
dim(FST_3Combined)

###Add a column that is identity to the AllmatchRADGTseq_FST_tags_to_keep
head(AllmatchRADGTseq_FST_tags_to_keep)
AllmatchRADGTseq_FST_tags_to_keep$Identity <- "NewRAD"
head(AllmatchRADGTseq_FST_tags_to_keep)

##########################
##### Create a data frame for each of the population groups within FST_3Combined data set
#with the loci sorted for the FST 

##Find the top FST for All
FST_3Combined_ALL <- FST_3Combined[,c("row.names","All", "Identity")]
head(FST_3Combined_ALL)
FST_3Combined_ALL <- FST_3Combined_ALL[order(FST_3Combined_ALL$All, decreasing=TRUE),]
head(FST_3Combined_ALL)

##Find the top FST for the Kusk
FST_3Combined_Kusk <- FST_3Combined[,c("row.names","Kusk", "Identity")]
head(FST_3Combined_Kusk)
FST_3Combined_Kusk <- FST_3Combined_Kusk[order(FST_3Combined_Kusk$Kusk, decreasing=TRUE),]
head(FST_3Combined_Kusk)

##Find the top FST for the Nush
FST_3Combined_Nush <- FST_3Combined[,c("row.names","Nush", "Identity")]
head(FST_3Combined_Nush)
FST_3Combined_Nush <- FST_3Combined_Nush[order(FST_3Combined_Nush$Nush, decreasing=TRUE),]
head(FST_3Combined_Nush)

##Find the top FST for the Kusk_Nush
FST_3Combined_Kusk_Nush <- FST_3Combined[,c("row.names","Kusk_Nush", "Identity")]
head(FST_3Combined_Kusk_Nush)
FST_3Combined_Kusk_Nush <- FST_3Combined_Kusk_Nush[order(FST_3Combined_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
head(FST_3Combined_Kusk_Nush)

########################################################################
##plot the range of the FST for each of the groups

AllmatchRADGTseq_FST_tags_to_keep_Ranked <- AllmatchRADGTseq_FST_tags_to_keep[order(AllmatchRADGTseq_FST_tags_to_keep$ALL_MRG, decreasing=TRUE),]
AllmatchRADGTseq_FST_tags_to_keep_Ranked$rank<-seq(1,dim(AllmatchRADGTseq_FST_tags_to_keep_Ranked)[1],by=1)
ggplot()+geom_point(data=AllmatchRADGTseq_FST_tags_to_keep_Ranked,aes(x=rank,y=ALL_MRG, color= Identity), size = .5)+ggtitle("FST range of newRAD ALL populations that match those in RADGTseq dataset ")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/AllmatchRADGTseq_FST_tags_to_keep_Ranked.pdf")


FST_3Combined_ALL_ranked <- FST_3Combined_ALL[order(FST_3Combined_ALL$All, decreasing=TRUE),]
FST_3Combined_ALL_ranked$rank<-seq(1,dim(FST_3Combined_ALL_ranked)[1],by=1)
ggplot()+geom_point(data=FST_3Combined_ALL_ranked,aes(x=rank,y=All, color= Identity), size = .5)+ggtitle("FST range of GTSeqONLY ALL populations ")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/AllmatchRADGTseq_FST_tags_to_keep_Ranked.pdf")

FST_3Combined_Kusk_ranked <- FST_3Combined_Kusk[order(FST_3Combined_Kusk$Kusk, decreasing=TRUE),]
FST_3Combined_Kusk_ranked$rank<-seq(1,dim(FST_3Combined_Kusk_ranked)[1],by=1)
ggplot()+geom_point(data=FST_3Combined_Kusk_ranked,aes(x=rank,y=Kusk,color = Identity), size = .5 )+ggtitle("FST range of GTSeqONLY Kuskokwim populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FST_3Combined_ALL_ranked.pdf")

FST_3Combined_Nush_ranked <- FST_3Combined_Nush[order(FST_3Combined_Nush$Nush, decreasing=TRUE),]
FST_3Combined_Nush_ranked$rank<-seq(1,dim(FST_3Combined_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=FST_3Combined_Nush_ranked,aes(x=rank,y=Nush ,color = Identity), size = .5 )+ggtitle("FST range of GTSeqONLY Nushagak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FST_3Combined_Nush_ranked.pdf")

FST_3Combined_Kusk_Nush_ranked <- FST_3Combined_Kusk_Nush[order(FST_3Combined_Kusk_Nush$Kusk_Nush, decreasing=TRUE),]
FST_3Combined_Kusk_Nush_ranked$rank<-seq(1,dim(FST_3Combined_Kusk_Nush_ranked)[1],by=1)
ggplot()+geom_point(data=FST_3Combined_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush,color = Identity), size = .5 )+ggtitle("FST range of GTSeqONLY Kuskokwim and Nushigak populations")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/ChinookPanel/FST_3Combined_Kusk_Nush_ranked.pdf")


#plot all five on the same graph in different colors
ggplot()+geom_point(data=FST_3Combined_ALL_ranked,aes(x=rank,y=All), color= "#93AA00") + #green
  geom_point(data=AllmatchRADGTseq_FST_tags_to_keep_Ranked,aes(x=rank,y=ALL_MRG), color= "magenta") + #magenta
  geom_point(data=FST_3Combined_Kusk_ranked,aes(x=rank,y=Kusk), color = "#00ADFA") + #blue
  geom_point(data=FST_3Combined_Nush_ranked,aes(x=rank,y=Nush), color= "#F8766D") + #salmon
  geom_point(data=FST_3Combined_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush), color = "#AE87FF") + #purple
  geom_hline(yintercept = 0, color = "black", size= .5) +  theme_bw() + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range of the five combined dataset population groups", x = "Rank", y = "FST", size = 20))  +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("#93AA00", "#E88527", "#00ADFA", "#F8766D","#AE87FF"))))

#plot all five  on the same graph, showing that the AllmatchRADGTseq_FST_tags_to_keep_Ranked is different than the FST_3Combined_ALL_ranked
ggplot()+geom_point(data=FST_3Combined_ALL_ranked,aes(x=rank,y=All, colour = "ALL")) +
  geom_point(data=AllmatchRADGTseq_FST_tags_to_keep_Ranked,aes(x=rank,y=ALL_MRG, color= "ALL_MRG")) +
  geom_point(data=FST_3Combined_Kusk_ranked,aes(x=rank,y=Kusk, colour = "Kusk")) +
  geom_point(data=FST_3Combined_Nush_ranked,aes(x=rank,y=Nush, colour = "Nush")) +
  geom_point(data=FST_3Combined_Kusk_Nush_ranked,aes(x=rank,y=Kusk_Nush, colour = "Kusk_Nush")) +
  scale_colour_manual(values = c("goldenrod1", "blue", "violetred1", "orangered1", "red1"), 
                      breaks = c("ALL","ALL_MRG", "Kusk", "Nush", "Kusk_Nush")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +  theme_bw() + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range of the four Combined dataset population groups", x = "Rank", y = "FST", size = 20))  +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("goldenrod1",  "blue", "violetred1", "orangered1", "red1"))))




##################################################
############## Select high FST loci  #############
##################################################

# FST_3Combined_ALL
# FST_3Combined_Kusk
# FST_3Combined_Nush
# FST_3Combined_Kusk_Nush
# AllmatchRADGTseq_FST_tags_to_keep_Ranked

################# Choose the top FST loci in each of the population groups for the data set 
#top FST loci for FST_3Combined_ALL_var
#change the variable to the number you want 
FST_3Combined_ALL_var <- 150 #<--------------------------------------------dictates the number of loci
FST_3Combined_ALL_set <- FST_3Combined_ALL[1:FST_3Combined_ALL_var,]
FST_3Combined_ALL_set_avgFST <- mean(FST_3Combined_ALL_set$All)
FST_3Combined_ALL_set_maxFST <- max(FST_3Combined_ALL_set$All)
FST_3Combined_ALL_set_minFST <- min(FST_3Combined_ALL_set$All)
cat("FST_3Combined_ALL FST max, min, average: ",FST_3Combined_ALL_set_maxFST, FST_3Combined_ALL_set_minFST, FST_3Combined_ALL_set_avgFST )

#top FST loci for FST_3Combined_Kusk
FST_3Combined_Kusk_var <- 200 #<--------------------------------------------dictates the number of loci
FST_3Combined_Kusk_set <- FST_3Combined_Kusk[1:FST_3Combined_Kusk_var,]
FST_3Combined_Kusk_set_avgFST <- mean(FST_3Combined_Kusk_set$Kusk)
FST_3Combined_Kusk_set_maxFST <- max(FST_3Combined_Kusk_set$Kusk)
FST_3Combined_Kusk_set_minFST <- min(FST_3Combined_Kusk_set$Kusk)
cat("FST_3Combined_Kusk FST max, min, average: ", FST_3Combined_Kusk_set_maxFST, FST_3Combined_Kusk_set_minFST, FST_3Combined_Kusk_set_avgFST )

#top FST loci for FST_3Combined_Nush
#change the variable to the number you want 
FST_3Combined_Nush_var <- 200 #<--------------------------------------------dictates the number of loci
FST_3Combined_Nush_set <- FST_3Combined_Nush[1:FST_3Combined_Nush_var,]
FST_3Combined_Nush_set_avgFST <- mean(FST_3Combined_Nush_set$Nush)
FST_3Combined_Nush_set_maxFST <- max(FST_3Combined_Nush_set$Nush)
FST_3Combined_Nush_set_minFST <- min(FST_3Combined_Nush_set$Nush)
cat("FST_3Combined_Nush FST max, min, average: ",FST_3Combined_Nush_set_maxFST, FST_3Combined_Nush_set_minFST, FST_3Combined_Nush_set_avgFST )

#top FST loci for FST_3Combined_Kusk_Nush
FST_3Combined_Kusk_Nush_var <- 275 #<--------------------------------------------dictates the number of loci
FST_3Combined_Kusk_Nush_set <- FST_3Combined_Kusk_Nush[1:FST_3Combined_Kusk_Nush_var,]
FST_3Combined_Kusk_Nush_set_avgFST <- mean(FST_3Combined_Kusk_Nush_set$Kusk_Nush)
FST_3Combined_Kusk_Nush_set_maxFST <- max(FST_3Combined_Kusk_Nush_set$Kusk_Nush)
FST_3Combined_Kusk_Nush_set_minFST <- min(FST_3Combined_Kusk_Nush_set$Kusk_Nush)
cat("FST_3Combined_Kusk_Nush FST max, min, average: ", FST_3Combined_Kusk_Nush_set_maxFST, FST_3Combined_Kusk_Nush_set_minFST, FST_3Combined_Kusk_Nush_set_avgFST )


#top FST loci for  AllmatchRADGTseq_FST_tags_to_keep_Ranked
AllmatchRADGTseq_FST_var <- 175 #<--------------------------------------------dictates the number of loci
AllmatchRADGTseq_FST_set <- AllmatchRADGTseq_FST_tags_to_keep_Ranked[1:AllmatchRADGTseq_FST_var,]
AllmatchRADGTseq_FST_set_avgFST <- mean(AllmatchRADGTseq_FST_set$ALL_MRG)
AllmatchRADGTseq_FST_set_maxFST <- max(AllmatchRADGTseq_FST_set$ALL_MRG)
AllmatchRADGTseq_FST_set_minFST <- min(AllmatchRADGTseq_FST_set$ALL_MRG)
cat("AllmatchRADGTseq FST max, min, average: ", AllmatchRADGTseq_FST_set_maxFST, AllmatchRADGTseq_FST_set_minFST, AllmatchRADGTseq_FST_set_avgFST )

### see which loci are in the five of the population groups
FST_3Combined_union <- c(FST_3Combined_Kusk_Nush_set$row.names, FST_3Combined_Nush_set$row.names, FST_3Combined_ALL_set$row.names, FST_3Combined_Kusk_set$row.names, AllmatchRADGTseq_FST_set$row.names)
FST_3Combined_union <- unique(FST_3Combined_union)
length(FST_3Combined_union)
head(FST_3Combined_union)

###Rank the loci 
## make a new dataframe that has the FST info for the union of the sets 
FST_3Combined_ALL_set$rank<-seq(1,dim(FST_3Combined_ALL_set)[1],by=1)
FST_3Combined_Kusk_set$rank<-seq(1,dim(FST_3Combined_Kusk_set)[1],by=1)
FST_3Combined_Nush_set$rank<-seq(1,dim(FST_3Combined_Nush_set)[1],by=1)
FST_3Combined_Kusk_Nush_set$rank<-seq(1,dim(FST_3Combined_Kusk_Nush_set)[1],by=1)
AllmatchRADGTseq_FST_set$rank<-seq(1,dim(AllmatchRADGTseq_FST_set)[1],by=1) 

#plot all four on the same graph in different colors
ggplot()+geom_point(data=FST_3Combined_ALL_set,aes(x=rank,y=All, colour = "ALL")) +
  geom_point(data=FST_3Combined_Kusk_set,aes(x=rank,y=Kusk, colour = "Kusk")) +
  geom_point(data=FST_3Combined_Nush_set,aes(x=rank,y=Nush, colour = "Nush")) +
  geom_point(data=FST_3Combined_Kusk_Nush_set,aes(x=rank,y=Kusk_Nush, colour = "Kusk_Nush")) +
  geom_point(data=AllmatchRADGTseq_FST_set,aes(x=rank,y=ALL_MRG, colour = "ALL_MRG")) +
  scale_colour_manual(values = c("goldenrod1", "violetred1", "orangered1", "red1", "blue"), 
                      breaks = c("ALL", "Kusk", "Nush", "Kusk_Nush", "ALL_MRG")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +  theme_bw() + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range of the four Combined dataset population groups", x = "Rank", y = "FST", size = 20))  +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("goldenrod1", "violetred1", "orangered1", "red1", "blue"))))

#plot all four on the same graph in different colors
ggplot()+geom_point(data=FST_3Combined_ALL_set,aes(x=rank,y=All, colour = "ALL")) +
  geom_point(data=FST_3Combined_Kusk_set,aes(x=rank,y=Kusk, colour = "Kusk")) +
  geom_point(data=FST_3Combined_Nush_set,aes(x=rank,y=Nush, colour = "Nush")) +
  geom_point(data=FST_3Combined_Kusk_Nush_set,aes(x=rank,y=Kusk_Nush, colour = "Kusk_Nush")) +
  geom_point(data=AllmatchRADGTseq_FST_set,aes(x=rank,y=ALL_MRG, colour = "ALL_MRG")) +
  scale_colour_manual(values = c("goldenrod1", "blue", "violetred1", "orangered1", "red1")) +
  geom_hline(yintercept = 0, color = "black", size= .5) +  theme_bw() + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range of the four Combined dataset population groups", x = "Rank", y = "FST", size = 20))  +
  guides(color = guide_legend(override.aes = list(shape = 19, size= 5, fill= c("goldenrod1", "violetred1", "orangered1", "red1", "blue"))))

#plot all four on the same graph in different colors
ggplot()+geom_point(data=FST_3Combined_ALL_set,aes(x=rank,y=All, colour = Identity), size = .5) +
  geom_point(data=FST_3Combined_Kusk_set,aes(x=rank,y=Kusk, colour = Identity), size = .5) +
  geom_point(data=FST_3Combined_Nush_set,aes(x=rank,y=Nush, colour = Identity), size = .5) +
  geom_point(data=FST_3Combined_Kusk_Nush_set,aes(x=rank,y=Kusk_Nush, colour = Identity), size = .5) +
  geom_point(data=AllmatchRADGTseq_FST_set,aes(x=rank,y=ALL_MRG, colour = Identity), size = .5) +
  geom_hline(yintercept = 0, color = "black", size= .5) +  theme_bw() + 
  theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "FST range of the four Combined dataset population groups", x = "Rank", y = "FST", size = 20))  

############
########Figure out how many of these are already in the panels- especially panel # 2, that tells the cook inlet pops apart 











########################################################
########################################################
########################################################
##I have not done the below ,but it is useful if we want to output some tables of these loci. 



### Write the FST table for the union of the GTseqONLY_FST_union_FST to a file 
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/GTseqONLY_FST_union.txt", "wb")
write.table(GTseqONLY_FST_union_FST,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)
## make a new dataframe that has the FST info for the union of the sets 
allCombined_FST_union_FST <- allCombined_FST_file[allCombined_FST_file$row.names %in% allCombined_FST_union,]  
allCombined_FST_union_FST_sort <-allCombined_FST_union_FST[order(allCombined_FST_union_FST$All, decreasing = TRUE),]
head(allCombined_FST_union_FST_sort)
dim(allCombined_FST_union_FST_sort)

### Write the FST table for the union of the allCombined_FST_union_FST to a file 
outputFile <- file("Z:/WORK/TARPEY/ChinookPanel/DataAnalysis/Genepop/allCombined_FST_union.txt", "wb")
write.table(allCombined_FST_union_FST,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

