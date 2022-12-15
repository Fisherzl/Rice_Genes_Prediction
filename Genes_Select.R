##############################################
#Author: Liu Zhu
#1.This part of the code is used to obtain the nodes shared by PPI network and co-expression network. 
#
#2.Building a seed node set of RWR-M.
#
#3.Building a test seed node set for LOOCV strategy.
##################################################
library(dplyr)
rm(list=ls())

###1###
#upload PPI network data.
ALL_PPI_data <- read.csv('ALL_PPI_Data_Thre_0.3.txt',sep = " ")
ALL_PPI_1 <- c(ALL_PPI_data[,1])
ALL_PPI_2 <- c(ALL_PPI_data[,2])
ALL_PPI_Proteins <- c(ALL_PPI_1,ALL_PPI_2)
Numbers_PPI_Proteins <- unique(ALL_PPI_Proteins)

#Constructing Co-expression network by different methods
#1. Constructing a gene co-expression network based on the Pearson correlation coefficient, 
#   with correlation-threshold > 0.5 as the screening criteria  
co_expression_data <- read.table("Differential_Relation_rnaseq0.5-co_expression.txt",sep = " ")

#2. Constructing a gene co-expression network based on the MIC(Maximal Information Coefficient), 
#   with correlation-threshold > 0.5 as the screening criteria
co_expression_data <- read.table("MIC_Relation_rnaseq0.5-co_expression.txt",sep = " ")

co_expression_data <- co_expression_data[,1:2]
co_expression_1 <- c(co_expression_data[,1])
co_expression_2 <- c(co_expression_data[,2])
ALL_co_Proteins <- c(co_expression_1,co_expression_2)
Numbers_co_Proteins <- unique(ALL_co_Proteins)

Intersection_Proteins <- c()

for (i in 1:length(Numbers_co_Proteins)) {
  if(Numbers_co_Proteins[i] %in% Numbers_PPI_Proteins){
    Intersection_Proteins <- append(Intersection_Proteins,Numbers_co_Proteins[i])
  }
}


#Screening proteins shared by ALL_PPI and co_expression
PPI_final_1 <- as.data.frame(matrix(nrow=0,ncol=0))
PPI_final_2 <- as.data.frame(matrix(nrow=0,ncol=0))
PPI_final_1 <- subset(ALL_PPI_data,ALL_PPI_data[,1] %in% Intersection_Proteins[,1])
PPI_final_2 <- subset(PPI_final_1,PPI_final_1[,2] %in% Intersection_Proteins[,1])
write.table(PPI_final_2,file = "PPI_2021_9_4.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)


co_final_1 <- as.data.frame(matrix(nrow=0,ncol=0))
co_final_2 <- as.data.frame(matrix(nrow=0,ncol=0))
co_final_1 <- subset(co_expression_data,co_expression_data[,1] %in% Intersection_Proteins[,1])
co_final_2 <- subset(co_final_1,co_final_1[,2] %in% Intersection_Proteins[,1])
write.table(co_final_2,file = "expression_2021_9_4.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)



#####################################################
#2.Building the seed node set of RWR-M.
####################################################

All_PPI <- read.csv('PPI_2021_9_4.txt',sep = " ",header = FALSE)
co_expression_data <- read.csv('expression_2021_9_4.txt',sep = " ",header = FALSE)

Yi_Protein_data <- read.csv('known_drought_stress_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')


Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% All_PPI[,1] | Yi_Protein_data[,1] %in% All_PPI[,2])
Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% co_expression_data[,1] | Select_Seed_Proteins_1[,1] %in% co_expression_data[,2])


write.table(Select_Seed_Proteins_2,file = "Intersection_Seed_Proteins.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)





#############################################################
#3.Building a test seed set for LOOCV strategy.
#############################################################
#Find seed nodes that are common to four different networks and use them as test sets to test the prediction effect.
MIC_Relation_expression_data <- read.csv("MIC_Relation_rnaseq0.6-co_expression.txt",sep = " ")
MIC_Relation_expression_data_1 <- c(MIC_Relation_expression_data[,1])
MIC_Relation_expression_data_2 <- c(MIC_Relation_expression_data[,2])
MIC_Relation_expression_data_Proteins <- c(MIC_Relation_expression_data_1,MIC_Relation_expression_data_2)
Numbers_MIC_Relation_expression_data_Proteins <- unique(MIC_Relation_expression_data_Proteins)

Diff_expression_data <- read.table("Differential_Relation_rnaseq0.5-co_expression.txt",sep = " ")
Diff_expression_data_1 <- c(Diff_expression_data[,1])
Diff_expression_data_2 <- c(Diff_expression_data[,2])
Diff_expression_data_Proteins <- c(Diff_expression_data_1,Diff_expression_data_2)
Numbers_Diff_expression_data_Proteins <- unique(Diff_expression_data_Proteins)


ALL_PPI_data <- read.table("DATA/ALL_PPI_Data.txt",sep = " ")
ALL_PPI_data_1 <- c(ALL_PPI_data[,1])
ALL_PPI_data_2 <- c(ALL_PPI_data[,2])
ALL_PPI_data_Proteins <- c(ALL_PPI_data_1,ALL_PPI_data_2)
Numbers_ALL_PPI_data_Proteins <- unique(ALL_PPI_data_Proteins)

STRING_data <- read.table("STRING_0.3.txt",sep = " ")
STRING_data_1 <- c(STRING_data[,1])
STRING_data_2 <- c(STRING_data[,2])
STRING_data_Proteins <- c(STRING_data_1,STRING_data_2)
Numbers_STRING_data_Proteins <- unique(STRING_data_Proteins)


Yi_Protein_data <- read.csv('known_drought_stress_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')


Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,1] | Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,2])
Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,1] | Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,2])
Select_Seed_Proteins_3 <- subset(Select_Seed_Proteins_2,Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,1] | Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,2])
Select_Seed_Proteins_4 <- subset(Select_Seed_Proteins_3,Select_Seed_Proteins_3[,1] %in% STRING_data[,1] | Select_Seed_Proteins_3[,1] %in% STRING_data[,2])

write.table(Select_Seed_Proteins_4,file = "Test_Seeds_brought_MIC0.3.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)


