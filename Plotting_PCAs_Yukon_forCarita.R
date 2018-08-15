### QC for the Yukon Project
### Individual Heterozygosity for the QC individuals and 
### Plotting the PCAS of the Yukon populations, QC and baseline
###   Checking to see where the QC cluster in relation to the populations in the baseline
### Carolyn Tarpey for Carita Pascal | August 2018
### ---------------------------------------

library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)




###############Import the pop info for the populations 
PolyGenCombined_POPMAP <- read.table("z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/191_newBaseline_polygen_combined_POPMAP.txt", header = TRUE, sep = '\t')
head(PolyGenCombined_POPMAP)


############Eigen tables #############
###Import the eigenvalues that were run in PLINK and merge them with the Pop info

newBaseline_polyGen_combined_table <- read.table("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/191_newBaseline_polyGen_combined_pca.eigenvec")
newBaseline_polyGen_combined_table <- merge(newBaseline_polyGen_combined_table, PolyGenCombined_POPMAP, by.x = 'V2', by.y = 'SampleName')
head(newBaseline_polyGen_combined_table)
unique(newBaseline_polyGen_combined_table$PopName)

newBaseline_polyGenQC_combined_table <- read.table("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/191_newBaseline_polyGenQC_combined_pca.eigenvec")
newBaseline_polyGenQC_combined_table <- merge(newBaseline_polyGenQC_combined_table, PolyGenCombined_POPMAP, by.x = 'V2', by.y = 'SampleName')
head(newBaseline_polyGenQC_combined_table)
unique(newBaseline_polyGenQC_combined_table$PopName)

newBaseline_subset_table <- read.table("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/191_newBaseline_subset_rel_pca.eigenvec")
newBaseline_subset_table <- merge(newBaseline_subset_table, PolyGenCombined_POPMAP, by.x = 'V2', by.y = 'SampleName')
head(newBaseline_subset_table)
unique(newBaseline_subset_table$PopName)

polyGenResults_haplotype_table <- read.table("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/191_polyGenResults_haplotype_pca.eigenvec")
polyGenResults_haplotype_table <- merge(polyGenResults_haplotype_table, PolyGenCombined_POPMAP, by.x = 'V2', by.y = 'SampleName')
head(polyGenResults_haplotype_table)
unique(polyGenResults_haplotype_table$SampleType)

polyGenResults_QCsubset_table <- read.table("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/191_polyGenResults_QCsubset_pca.eigenvec")
polyGenResults_QCsubset_table <- merge(polyGenResults_QCsubset_table, PolyGenCombined_POPMAP, by.x = 'V2', by.y = 'SampleName')
head(polyGenResults_QCsubset_table)
unique(polyGenResults_QCsubset_table$SampleType)

NA_col_4<-  c("KANVI" = "magenta", "KSALC"="#7297ee","KWHITERA" = "#00158a", "KLYUTF"="#e0957e")
All_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Susitna" = "#8FA9B7", 
            "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_3<-  c("NORM" = "#AAAAD4", "OQC"="#5555AA","QC" = "cyan")
NA_col_2<-  c( "OQC"="#5555AA","QC" = "cyan")

################# PCs one and two ###############################
pdf("Z:/WORK/TARPEY/ExamplesForOthers/YukonPCA/PCAs/newBaseline_polyGen_combined.pdf", width = 9, height = 7)

ggplot(data = newBaseline_polyGen_combined_table) + geom_point(aes(x =V3, y = V4,  color = PopName, shape = SampleType), alpha = .8, size = 2) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "newBaseline_polyGen PC #s 1 & 2", x = "newBaseline_polyGen 1", y = "newBaseline_polyGen 2", size = 20))+ 
  scale_color_manual(breaks = c("KANVI","KSALC","KWHITERA","KLYUTF"), values = NA_col_4) 

ggplot(data = newBaseline_polyGenQC_combined_table) + geom_point(aes(x =V3, y = V4,  color = PopName, shape = SampleType), alpha = .8, size = 2) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "newBaseline_polyGenQC_combined PC #s 1 & 2", x = "newBaseline_polyGenQC_combined 1", y = "newBaseline_polyGenQC_combined 2 ", size = 20))+ 
  scale_color_manual(breaks = c("KANVI","KSALC","KWHITERA","KLYUTF"), values = NA_col_4) 


ggplot(data = newBaseline_subset_table) + geom_point(aes(x =V3, y = V4,  color = PopName, shape = SampleType), alpha = .8, size = 2) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "newBaseline_subset PC #s 1 & 2", x = "newBaseline_subset 1", y = "newBaseline_subset 2", size = 20))+ 
  scale_color_manual(breaks = c("KANVI","KSALC","KWHITERA","KLYUTF"), values = NA_col_4) 

ggplot(data = polyGenResults_haplotype_table) + geom_point(aes(x =V3, y = V4,  color = SampleType), alpha = .8, size = 2) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "polyGenResults_haplotype PC #s 1 & 2", x = "polyGenResults_haplotype 1", y = "polyGenResults_haplotype 2", size = 20))+ 
  scale_color_manual(breaks = c("NORM","OQC","QC"), values = NA_col_3) 

ggplot(data = polyGenResults_QCsubset_table) + geom_point(aes(x =V3, y = V4,  color = SampleType, shape = SampleType), alpha = .8, size = 2) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "polyGenResults_QCsubset PC #s 1 & 2", x = "polyGenResults_QCsubset 1", y = "polyGenResults_QCsubset 2", size = 20))+ 
  scale_color_manual(breaks = c("OQC","QC"), values = NA_col_2) 

dev.off()





