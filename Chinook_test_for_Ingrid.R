###  Chinook panel: editing the raw genepop file as a test for Ingrid
###    
###    
### Carolyn Tarpey | July 2018 
### ---------------------------------------

### Read in the data files we will need ########

#load the raw genepop file, edited for R. Includes all individuals and loci:  
chin_raw_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/batch_13.allCombined_NoCookInlet_R_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_raw_genepop)

chin_raw_genepop_edit <- chin_raw_genepop[1:200,1:200]
dim(chin_raw_genepop_edit)
head(chin_raw_genepop_edit)

colnames(chin_raw_genepop_edit)<-gsub("X","",colnames(chin_raw_genepop_edit))


write.table(chin_raw_genepop_edit,"Z:/WORK/TARPEY/ExamplesForOthers/batch_13.allCombined_200SNPS_200Inds_genepop_R.txt",quote=FALSE,row.names=TRUE, sep= "\t")
