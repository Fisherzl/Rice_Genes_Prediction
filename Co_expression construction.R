#########################################################
#Author: Liu Zhu
#This part is used to build a gene co-expression network,
#including MIC and Pearson correlation coefficient.
########################################################
rm(list = ls())
library(Hmisc)
library(WGCNA)
library(minerva)

###############################################
#1¡¢Starting from the expression matrix after differential expression screening, 
#this data analyzes the co-expression network established under the 0.5, 0.6 and 0.7 thresholds.

#DEGs
mydata <- read.csv("Rna-seq_Differential_expression_Select.csv")
rownames(mydata) <- mydata$X
mydata <- mydata[,-(1)]
mydata <- t(mydata)
res <- cor(mydata)
round(res, 2)

probes <- colnames(res)
dimnames(res) <- list(probes,probes)
cyt <- exportNetworkToCytoscape(res,
                                edgeFile = paste("exprMAt",".edge.txt",sep = ""),
                                nodeFile = paste("exprMAt",".node.txt",sep = ""),
                                weighted = TRUE,threshold = 0,
                                nodeNames = probes)
woyaode <- cyt[["edgeData"]]
woyaode <- woyaode[,1:3]
woyaode <- subset(woyaode,woyaode$weight>0.5)

write.table(woyaode,file = "Differential_Relation_rnaseq0.5-co_expression.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)



#2¡¢30% of the genes with high phenotypic correlation were screened by MIC, 
#and the correlation between genes was discovered by MIC to construct a gene co-expression network

mydata <- read.csv("myrice_drought_rnamic30.csv")
rownames(mydata) <- mydata$X
mydata <- mydata[,-(1:2)]
MIC30_matrix <- mine(mydata)

MIC30_matrix_m <- as.matrix(MIC30_matrix)

MIC30_matrix_M <- matrix(0,nrow = 29521,ncol = 1)
MIC30_matrix_M <- MIC30_matrix_m[1,]
MIC30_matrix_M_MIC <- as.matrix(MIC30_matrix_M[["MIC"]])
rownames(MIC30_matrix_M_MIC) <- colnames(mydata)



probes <- colnames(MIC30_matrix_M_MIC[,-1])
dimnames(MIC30_matrix_M_MIC[,-1]) <- list(probes,probes)
cyt <- exportNetworkToCytoscape(MIC30_matrix_M_MIC[,-1],
                                edgeFile = paste("exprMAt",".edge.txt",sep = ""),
                                nodeFile = paste("exprMAt",".node.txt",sep = ""),
                                weighted = TRUE,threshold = 0,
                                nodeNames = probes)
woyaode <- cyt[["edgeData"]]
woyaode <- woyaode[,1:3]
woyaode <- subset(woyaode,woyaode$weight > 0.5)

write.table(woyaode,file = "F:/xuejie_data/0.5_Relation_rnaseq_Differexpression_MIC_Select.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)





