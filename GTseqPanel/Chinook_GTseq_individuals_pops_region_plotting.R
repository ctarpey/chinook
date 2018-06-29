###  Chinook panel: Individual and Pop demographic plotting
###    Exploratory analysis to visualize the locations and population numbers
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


PopInfo_RAD_data <- read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/populationMap_NoCookInlet.txt", header =TRUE)
head(PopInfo_RAD_data)
dim(PopInfo_RAD_data)


PopInfo_GT_data <- read.delim("Z:/WORK/TARPEY/ChinookPanel/combined_RAD_GTseq_data/PopulationMap_847RAD_taq_genepop.txt", header =TRUE)
head(PopInfo_GT_data)
dim(PopInfo_GT_data)

PopInfo_GTONLY_data <- read.delim("Z:/WORK/TARPEY/ChinookPanel/GTseqOnlyData/PopulationMap_GTseqDataONLY.txt", header =TRUE)
head(PopInfo_GTONLY_data)
dim(PopInfo_GTONLY_data)

######################### Population Numbers for the NoCookInlet genepop file, pre-filtering loci and individuals###############

#Plot the difference in the population group numbers 

#summary of all the individuals by population
RAD_Inds_by_pop<- table(PopInfo_RAD_data[,2])
RAD_Inds_by_pop <- as.data.frame(RAD_Inds_by_pop)
colnames(RAD_Inds_by_pop) <- c("Population","NumInds")
head(RAD_Inds_by_pop)
dim(RAD_Inds_by_pop)

#Number of individuals in each of the scales of reporting groups
RAD_VeryBroad <- table(PopInfo_RAD_data[,3])
RAD_VeryBroad <- as.data.frame(RAD_VeryBroad)
colnames(RAD_VeryBroad) <- c("ReportingGroup","NumInds")
head(RAD_VeryBroad)
dim(RAD_VeryBroad)

RAD_Broad <- table(PopInfo_RAD_data[,5])
RAD_Broad <- as.data.frame(RAD_Broad)
colnames(RAD_Broad) <- c("ReportingGroup","NumInds")
head(RAD_Broad)
dim(RAD_Broad)

RAD_Fine <- table(PopInfo_RAD_data[,6])
RAD_Fine <- as.data.frame(RAD_Fine)
colnames(RAD_Fine) <- c("ReportingGroup","NumInds")
head(RAD_Fine)
dim(RAD_Fine)

#################Colors for the NoCookInlet genepop file, pre-filtering loci and individuals ##########
Chin_Pops_16<- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#5aae61",'#a6cee3',"#3288bd",
                    '#cab2d6','#5e4fa2',"#762a83",'#2d004b')
plot(rep(1,16),col=Chin_Pops_16,pch=19,cex=7)

Chin_Pops_7<- c("#d53e4f","#f46d43","#fdae61","#abdda4","#66c2a5","#3288bd",'#5e4fa2')
plot(rep(1,7),col=Chin_Pops_7,pch=19,cex=7)

Chin_Pops_4<- c("#d53e4f","#fdae61","#3288bd",'#5e4fa2')
plot(rep(1,4),col=Chin_Pops_4,pch=19,cex=7)

Chin_Pops_2<- c("#fdae61","#3288bd")
plot(rep(1,2),col=Chin_Pops_2,pch=19,cex=17)

# colfunc <- colorRampPalette(c("#66c2a5", "#5e4fa2"))
# susitna_3 <- colfunc(3)
# plot(rep(1,3),col=susitna_3,pch=19,cex=7)


#################Plots for the NoCookInlet genepop file, pre-filtering loci and individuals####################

#plot of all the individuals in each population 
ggplot(data=RAD_Inds_by_pop, aes(RAD_Inds_by_pop$Population, RAD_Inds_by_pop$NumInds)) + theme_bw() +
  geom_bar(aes(fill= RAD_Inds_by_pop$Population),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Population") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_16) +
  geom_text(aes(x = RAD_Inds_by_pop$Population, y = RAD_Inds_by_pop$NumInds, label = RAD_Inds_by_pop$NumInds, 
  group = RAD_Inds_by_pop$Population), position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + 
  scale_y_continuous(limits = c(0, 80))

#plot of all the individuals in each very broad group
ggplot(data=RAD_VeryBroad, aes(RAD_VeryBroad$ReportingGroup, RAD_VeryBroad$NumInds)) + theme_bw() +
  geom_bar(aes(fill= RAD_VeryBroad$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Very Broad Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Very Broad Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_2)  + 
  geom_text(aes(x = RAD_VeryBroad$ReportingGroup, y = RAD_VeryBroad$NumInds, label = RAD_VeryBroad$NumInds, group = RAD_VeryBroad$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 700))

#plot of all the individuals in each broad group
ggplot(data=RAD_Broad, aes(RAD_Broad$ReportingGroup, RAD_Broad$NumInds)) + theme_bw() +
  geom_bar(aes(fill= RAD_Broad$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Broad Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Broad Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_4)  + 
  geom_text(aes(x = RAD_Broad$ReportingGroup, y = RAD_Broad$NumInds, label = RAD_Broad$NumInds, group = RAD_Broad$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 450))


#plot of all the individuals in each fine reporting group
ggplot(data=RAD_Fine, aes(RAD_Fine$ReportingGroup, RAD_Fine$NumInds)) + theme_bw() +
  geom_bar(aes(fill= RAD_Fine$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Fine Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Fine Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_7)  + 
  geom_text(aes(x = RAD_Fine$ReportingGroup, y = RAD_Fine$NumInds, label = RAD_Fine$NumInds, group = RAD_Fine$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 400))




######################### Population Numbers for the 847 RAD/GTseq file, pre-filtering ###############

#Plot the difference in the population group numbers 

#summary of all the individuals by population
GTseq_Inds_by_pop<- table(PopInfo_GT_data[,2])
GTseq_Inds_by_pop <- as.data.frame(GTseq_Inds_by_pop)
colnames(GTseq_Inds_by_pop) <- c("Population","NumInds")
head(GTseq_Inds_by_pop)
dim(GTseq_Inds_by_pop)


GTseq_Broad <- table(PopInfo_GT_data[,4])
GTseq_Broad <- as.data.frame(GTseq_Broad)
colnames(GTseq_Broad) <- c("ReportingGroup","NumInds")
head(GTseq_Broad)
dim(GTseq_Broad)

GTseq_Fine <- table(PopInfo_GT_data[,5])
GTseq_Fine <- as.data.frame(GTseq_Fine)
colnames(GTseq_Fine) <- c("ReportingGroup","NumInds")
head(GTseq_Fine)
dim(GTseq_Fine)

#################Colors for the 847 RAD/GTseq file, pre-filtering ##########
Chin_Pops_18<- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#5aae61",'#a6cee3',"#3288bd",
                 "#253494",'#cab2d6','#5e4fa2',"#762a83","#ae017e",'#2d004b')
plot(rep(1,18),col=Chin_Pops_18,pch=19,cex=7)

Chin_Pops_3<- c("#3288bd","#fdae61",'#5e4fa2')
plot(rep(1,3),col=Chin_Pops_3,pch=19,cex=17)

Chin_Pops_4<- c("#d53e4f","#3288bd","#fdae61",'#5e4fa2')
plot(rep(1,4),col=Chin_Pops_4,pch=19,cex=7)

# colfunc <- colorRampPalette(c("#66c2a5", "#5e4fa2"))
# susitna_3 <- colfunc(3)
# plot(rep(1,3),col=susitna_3,pch=19,cex=7)


#################Plots for the 847 RAD/GTseq file, pre-filtering####################

#plot of all the individuals in each population 
ggplot(data=GTseq_Inds_by_pop, aes(GTseq_Inds_by_pop$Population, GTseq_Inds_by_pop$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseq_Inds_by_pop$Population),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Population") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_18) +
  geom_text(aes(x = GTseq_Inds_by_pop$Population, y = GTseq_Inds_by_pop$NumInds, label = GTseq_Inds_by_pop$NumInds, 
                group = GTseq_Inds_by_pop$Population), position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + 
  scale_y_continuous(limits = c(0, 110))

#plot of all the individuals in each very broad group
ggplot(data=GTseq_Broad, aes(GTseq_Broad$ReportingGroup, GTseq_Broad$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseq_Broad$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Broad Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Broad Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_3)  + 
  geom_text(aes(x = GTseq_Broad$ReportingGroup, y = GTseq_Broad$NumInds, label = GTseq_Broad$NumInds, group = GTseq_Broad$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 900))

#plot of all the individuals in each fine reporting group
ggplot(data=GTseq_Fine, aes(GTseq_Fine$ReportingGroup, GTseq_Fine$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseq_Fine$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Fine Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Fine Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_4)  + 
  geom_text(aes(x = GTseq_Fine$ReportingGroup, y = GTseq_Fine$NumInds, label = GTseq_Fine$NumInds, group = GTseq_Fine$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 800))



######################### Population Numbers for the 1092 GTseq ONLY file, pre-filtering ###############

#Plot the difference in the population group numbers 

#summary of all the individuals by population
GTseqONLY_Inds_by_pop<- table(PopInfo_GTONLY_data[,2])
GTseqONLY_Inds_by_pop <- as.data.frame(GTseqONLY_Inds_by_pop)
colnames(GTseqONLY_Inds_by_pop) <- c("Population","NumInds")
head(GTseqONLY_Inds_by_pop)
dim(GTseqONLY_Inds_by_pop)

GTseqONLY_Broad <- table(PopInfo_GTONLY_data[,5])
GTseqONLY_Broad <- as.data.frame(GTseqONLY_Broad)
colnames(GTseqONLY_Broad) <- c("ReportingGroup","NumInds")
head(GTseqONLY_Broad)
dim(GTseqONLY_Broad)

GTseqONLY_Fine <- table(PopInfo_GTONLY_data[,6])
GTseqONLY_Fine <- as.data.frame(GTseqONLY_Fine)
colnames(GTseqONLY_Fine) <- c("ReportingGroup","NumInds")
head(GTseqONLY_Fine)
dim(GTseqONLY_Fine)

#################Colors for the 1092 GTseq ONLY file, pre-filtering ##########
Chin_Pops_17<- c("#9e0142","#d53e4f","#d53e4f","#f46d43","#fdae61","#fee08b","#fee08b","#e6f598","#abdda4","#66c2a5","#5aae61",'#a6cee3',"#3288bd",
                 "#253494",'#cab2d6','#5e4fa2',"#762a83","#ae017e",'#2d004b')
plot(rep(1,19),col=Chin_Pops_17,pch=19,cex=7)


Chin_Pops_3<- c("#3288bd","#fdae61",'#5e4fa2')
plot(rep(1,3),col=Chin_Pops_3,pch=19,cex=17)

Chin_Pops_4<- c("#d53e4f","#3288bd","#fdae61",'#5e4fa2')
plot(rep(1,4),col=Chin_Pops_4,pch=19,cex=7)

# colfunc <- colorRampPalette(c("#66c2a5", "#5e4fa2"))
# susitna_3 <- colfunc(3)
# plot(rep(1,3),col=susitna_3,pch=19,cex=7)


#################Plots for the 1092 GTseq ONLY file, pre-filtering####################

#plot of all the individuals in each population 
ggplot(data=GTseqONLY_Inds_by_pop, aes(GTseqONLY_Inds_by_pop$Population, GTseqONLY_Inds_by_pop$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseqONLY_Inds_by_pop$Population),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Population") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Population") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_17) +
  geom_text(aes(x = GTseqONLY_Inds_by_pop$Population, y = GTseqONLY_Inds_by_pop$NumInds, label = GTseqONLY_Inds_by_pop$NumInds, 
                group = GTseqONLY_Inds_by_pop$Population), position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + 
  scale_y_continuous(limits = c(0, 110))

#plot of all the individuals in each very broad group
ggplot(data=GTseqONLY_Broad, aes(GTseqONLY_Broad$ReportingGroup, GTseqONLY_Broad$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseqONLY_Broad$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Broad Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Broad Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_3)  + 
  geom_text(aes(x = GTseqONLY_Broad$ReportingGroup, y = GTseqONLY_Broad$NumInds, label = GTseqONLY_Broad$NumInds, group = GTseqONLY_Broad$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 550))

#plot of all the individuals in each fine reporting group
ggplot(data=GTseqONLY_Fine, aes(GTseqONLY_Fine$ReportingGroup, GTseqONLY_Fine$NumInds)) + theme_bw() +
  geom_bar(aes(fill= GTseqONLY_Fine$ReportingGroup),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Fine Reporting Group") + ylab("Number of Individuals") + ggtitle("Number of Individuals per Fine Reporting Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Chin_Pops_4)  + 
  geom_text(aes(x = GTseqONLY_Fine$ReportingGroup, y = GTseqONLY_Fine$NumInds, label = GTseqONLY_Fine$NumInds, group = GTseqONLY_Fine$ReportingGroup),
            position = position_dodge(width = .9), vjust = -0.5, size = 3.5) + scale_y_continuous(limits = c(0, 500))





############################# LOCUS STATS############################
#These are for when there is a filtered set that has the statistics 


################################## Statistics Per Data Set

### Calculate Individual Genotype Rate per Population 
dist <-dim(Genepop)[2]
Genepop_IGR_temp<-apply(Genepop[,10:dist],1,function(x) 1-(sum(x=="0000")/n_loci))
Genepop_IGR<-cbind(Genepop_IGR_temp,Genepop, by = "Names")
colnames(Genepop_IGR)[1]<- "SampleGenoRate"
dim(Genepop_IGR)
Genepop_IGR[1:10,1:10]

######## Calculate Locus Genotype Rate per Population 
#create a dataframe for our Locus Geno Rate results
Genepop_LGR <- matrix(nrow=n_pops, ncol=n_loci) 
rownames(Genepop_LGR)<-as.vector(unique(Genepop$TRUENAME)) #name the rows by the TRUE population names
colnames(Genepop_LGR)<-loci_names
dim(Genepop_LGR)

QC_test <-vector()

for (i in 1:n_pops){
  tempPop<-(unique(Genepop$TRUENAME)[i])
  #QC_test[i]<- tempPop
  tempGeno<-subset(Genepop, TRUENAME == tempPop)
  tempLGR<-apply(tempGeno[,10:dim(Genepop)[2]],2,function(x) 1-(sum(x=="0000")/dim(tempGeno)[1]))
  Genepop_LGR[i,]<-tempLGR
}

Genepop_LGR[1:15,1:15]
dim(Genepop_LGR)


################# Calculate MAF per population 

#get minor allele frequency for loci
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

#create a dataframe for our MAF results

Genepop_MAF<-matrix(nrow=n_pops, ncol=n_loci) 
rownames(Genepop_MAF)<-as.vector(unique(Genepop$TRUENAME)) #name the rows by the population names
colnames(Genepop_MAF)<-loci_names

for (i in 1:n_pops){
  tempPop<-(unique(Genepop$TRUENAME)[i])
  tempGeno<-subset(Genepop, TRUENAME == tempPop)
  tempMAF<-apply(tempGeno[,10:dim(Genepop)[2]],2,calculateMAF)
  Genepop_MAF[i,]<-tempMAF
}

Genepop_MAF[1:5,1:5]
##################### Run Genepop with the genepop file for HWE and Fis and import the results here ////// This was taken from SecondPinkFiltering.R

############ Heterozygosity per Population
## This is only run on the Completely filtered 23,759 loci and it was calculated with a perl script written by garrett: countHets_genepop.pl

HetCounts <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/Genepop/PercentHetsPerInd.txt" ,
                        stringsAsFactors = FALSE,header= TRUE)
head(HetCounts)
dim(HetCounts)

#add the PopINFO to the het counts
HetCounts<-merge(HetCounts,PED, by= "Names")
dim(HetCounts)
HetCounts[1:15,1:5]
HetCounts<-merge(HetCounts,POP_INFO, by= "POP")
dim(HetCounts)
HetCounts[1:15,1:12]

############################################## PLOTTING The statistics per population  
xx <- ggplot(a) + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Show the colors used for a default palette of 18
show_col(hue_pal()(18))

# this code is from Eleni, she used it to plot something similar in her data. 
#This code requires a dataframe with the individuals and a column that has population designation
#ggplot(mydata, aes(x=loci, fill= population)) + geom_histogram(data= mydata, bins = 20) + facet_wrap(~population) 

####################### Individual Genotype Rate
Genepop_IGR[1:5,1:15]


# plot Individual Genotype Rate Post_Filtered_data
ggplot(Genepop_IGR, aes(x=SampleGenoRate, fill= TRUENAME)) + 
  geom_histogram(data= Genepop_IGR, bins = 20) + facet_wrap(~TRUENAME) + theme_bw()  + guides(fill= FALSE) +
  ggtitle("Counts of Genotype Rate per Individual by Population")


####################### Locus Genotype Rate
Genepop_LGR[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(Genepop_LGR["AMUR10",] )
b <-as.data.frame(Genepop_LGR["LAKEL07",] )
c <-as.data.frame(Genepop_LGR["SUSIT14",] )
d <-as.data.frame(Genepop_LGR["NOME91",] )
e <-as.data.frame(Genepop_LGR["NOME94",] )
f <-as.data.frame(Genepop_LGR["SNOH03",] )
g <-as.data.frame(Genepop_LGR["SNOH96",] )
h <-as.data.frame(Genepop_LGR["LAKEL06",] )
i <-as.data.frame(Genepop_LGR["TAUY09",] )
j <-as.data.frame(Genepop_LGR["TAUY12",] )
k <-as.data.frame(Genepop_LGR["AMUR11",] )
l <-as.data.frame(Genepop_LGR["SUSIT13",] )
m <-as.data.frame(Genepop_LGR["HAYLY09",] )
n <-as.data.frame(Genepop_LGR["HAYLY10",] )
o <-as.data.frame(Genepop_LGR["KOPPE91",] )
p <-as.data.frame(Genepop_LGR["KOPPE96",] )
q <-as.data.frame(Genepop_LGR["KUSHI06",] )
r <-as.data.frame(Genepop_LGR["KUSHI07",] )

#get the population names in order
rownames(Genepop_LGR)

#assign each of the populations a ggplot to be called later
a1<- xx + geom_histogram(data= a, aes(x=a), bins = 50, fill="#F8766D") +  ggtitle("AMUR10") 
b1<- xx + geom_histogram(data= b, aes(x=b), bins = 50, fill="#E88526") +  ggtitle("LAKEL07") 
c1<- xx + geom_histogram(data= c, aes(x=c), bins = 50, fill="#D39200") +  ggtitle("SUSIT14") 
d1<- xx + geom_histogram(data= d, aes(x=d), bins = 50, fill="#B79F00") +  ggtitle("NOME91") 
e1<- xx + geom_histogram(data= e, aes(x=e), bins = 50, fill="#93AA00") +  ggtitle("NOME94") 
f1<- xx + geom_histogram(data= f, aes(x=f), bins = 50, fill="#5EB300") +  ggtitle("SNOH03") 
g1<- xx + geom_histogram(data= g, aes(x=g), bins = 50, fill="#00BA38") +  ggtitle("SNOH96") 
h1<- xx + geom_histogram(data= h, aes(x=h), bins = 50, fill="#00BF74") +  ggtitle("LAKEL06") 
i1<- xx + geom_histogram(data= i, aes(x=i), bins = 50, fill="#00C19F") +  ggtitle("TAUY09") 
j1<- xx + geom_histogram(data= j, aes(x=j), bins = 50, fill="#00BFC4") +  ggtitle("TAUY12") 
k1<- xx + geom_histogram(data= k, aes(x=k), bins = 50, fill="#00B9C3") +  ggtitle("AMUR11") 
l1<- xx + geom_histogram(data= l, aes(x=l), bins = 50, fill="#00ADFA") +  ggtitle("SUSIT13") 
m1<- xx + geom_histogram(data= m, aes(x=m), bins = 50, fill="#619CFF") +  ggtitle("HAYLY09") 
n1<- xx + geom_histogram(data= n, aes(x=n), bins = 50, fill="#AE87FF") +  ggtitle("HAYLY10") 
o1<- xx + geom_histogram(data= o, aes(x=o), bins = 50, fill="#DB72FB") +  ggtitle("KOPPE91") 
p1<- xx + geom_histogram(data= p, aes(x=p), bins = 50, fill="#F564E3") +  ggtitle("KOPPE96") 
q1<- xx + geom_histogram(data= q, aes(x=q), bins = 50, fill="#FF61C3") +  ggtitle("KUSHI06") 
r1<- xx + geom_histogram(data= r, aes(x=r), bins = 50, fill="#FF699C") +  ggtitle("KUSHI07") 

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Counts of Locus Genotype Rate by Population")

####################### Minor Allele Frequency BY pop
Genepop_MAF[1:5,1:15]

#Subset the matrix to get each population's numbers
a <-as.data.frame(Genepop_MAF["AMUR10",] )
b <-as.data.frame(Genepop_MAF["LAKEL07",] )
c <-as.data.frame(Genepop_MAF["SUSIT14",] )
d <-as.data.frame(Genepop_MAF["NOME91",] )
e <-as.data.frame(Genepop_MAF["NOME94",] )
f <-as.data.frame(Genepop_MAF["SNOH03",] )
g <-as.data.frame(Genepop_MAF["SNOH96",] )
h <-as.data.frame(Genepop_MAF["LAKEL06",] )
i <-as.data.frame(Genepop_MAF["TAUY09",] )
j <-as.data.frame(Genepop_MAF["TAUY12",] )
k <-as.data.frame(Genepop_MAF["AMUR11",] )
l <-as.data.frame(Genepop_MAF["SUSIT13",] )
m <-as.data.frame(Genepop_MAF["HAYLY09",] )
n <-as.data.frame(Genepop_MAF["HAYLY10",] )
o <-as.data.frame(Genepop_MAF["KOPPE91",] )
p <-as.data.frame(Genepop_MAF["KOPPE96",] )
q <-as.data.frame(Genepop_MAF["KUSHI06",] )
r <-as.data.frame(Genepop_MAF["KUSHI07",] )

#get the population names in order
rownames(Genepop_MAF)

#assign each of the populations a ggplot to be called later
a1<- xx + geom_histogram(data= a, aes(x=a), bins = 50, fill="#F8766D") +  ggtitle("AMUR10") 
b1<- xx + geom_histogram(data= b, aes(x=b), bins = 50, fill="#E88526") +  ggtitle("LAKEL07") 
c1<- xx + geom_histogram(data= c, aes(x=c), bins = 50, fill="#D39200") +  ggtitle("SUSIT14") 
d1<- xx + geom_histogram(data= d, aes(x=d), bins = 50, fill="#B79F00") +  ggtitle("NOME91") 
e1<- xx + geom_histogram(data= e, aes(x=e), bins = 50, fill="#93AA00") +  ggtitle("NOME94") 
f1<- xx + geom_histogram(data= f, aes(x=f), bins = 50, fill="#5EB300") +  ggtitle("SNOH03") 
g1<- xx + geom_histogram(data= g, aes(x=g), bins = 50, fill="#00BA38") +  ggtitle("SNOH96") 
h1<- xx + geom_histogram(data= h, aes(x=h), bins = 50, fill="#00BF74") +  ggtitle("LAKEL06") 
i1<- xx + geom_histogram(data= i, aes(x=i), bins = 50, fill="#00C19F") +  ggtitle("TAUY09") 
j1<- xx + geom_histogram(data= j, aes(x=j), bins = 50, fill="#00BFC4") +  ggtitle("TAUY12") 
k1<- xx + geom_histogram(data= k, aes(x=k), bins = 50, fill="#00B9C3") +  ggtitle("AMUR11") 
l1<- xx + geom_histogram(data= l, aes(x=l), bins = 50, fill="#00ADFA") +  ggtitle("SUSIT13") 
m1<- xx + geom_histogram(data= m, aes(x=m), bins = 50, fill="#619CFF") +  ggtitle("HAYLY09") 
n1<- xx + geom_histogram(data= n, aes(x=n), bins = 50, fill="#AE87FF") +  ggtitle("HAYLY10") 
o1<- xx + geom_histogram(data= o, aes(x=o), bins = 50, fill="#DB72FB") +  ggtitle("KOPPE91") 
p1<- xx + geom_histogram(data= p, aes(x=p), bins = 50, fill="#F564E3") +  ggtitle("KOPPE96") 
q1<- xx + geom_histogram(data= q, aes(x=q), bins = 50, fill="#FF61C3") +  ggtitle("KUSHI06") 
r1<- xx + geom_histogram(data= r, aes(x=r), bins = 50, fill="#FF699C") +  ggtitle("KUSHI07") 

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1, nrow=4, top="Counts of Minor Allele Frequency by Population")


####################### Individual Heterozygosity by Population
HetCounts[1:5,1:5]

# plot Counts of Individual Heterozygosity by Population
ggplot(HetCounts, aes(x=PercentHeterozygous, fill= TRUENAME)) + 
  geom_histogram(data=HetCounts, bins = 20) + facet_wrap(~TRUENAME) + theme_bw()  + guides(fill= FALSE) +
  ggtitle("Counts of Individual Percent Heterozygosity by Population")

#Mean Individual Heterozygosity by Population

# plot Mean Individual Heterozygosity by Population
ggplot(HetCounts, aes(x=factor(TRUENAME), y=PercentHeterozygous, fill=TRUENAME)) + stat_summary(fun.y= "mean", geom="bar") + theme_bw() + 
  ggtitle("Mean Individual Percent Heterozygosity by Population") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)


######################################

#####THERE ARE MORE INDIVIDUAL POPULATION AND GROUP SAMPLE NUMBER PLOTS IN R CODE Pink_PLINK_LD_popPlots.R


####################### Number of individuals by population PRE AND POST FILTERING 
popmap465_NEW<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap465_NEWNAMES_NEWSAMPLES.txt")
popmap492<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap_NEWNAMES.txt")

head(popmap465_NEW)
head(popmap492)

colnames(popmap465_NEW) <- c("Sample","Pop")
colnames(popmap492) <- c("Sample","Pop")

# Post Filtering
ggplot(data=popmap465_NEW) + geom_bar(aes(x=popmap465_NEW$Pop, fill=popmap465_NEW$Pop), stat= "count") + theme_bw() + 
  ggtitle("Number of Samples per Population, After Filtering") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)

#Pre-Filtering
ggplot(data=popmap492) + geom_bar(aes(x=popmap492$Pop, fill=popmap492$Pop), stat= "count") + theme_bw() + 
  ggtitle("Number of Samples per Population, Before Filtering") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)




