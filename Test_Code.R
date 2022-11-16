
#处理共表达网络和PPI网络的共同蛋白质
rm(list=ls())
###1###
ALL_PPI_data <- read.csv('E:/zhuliu_File/R Code_PPI/DATA/1-Threshold Proteins/ALL_PPI_Data_Thre_0.3.txt',sep = " ")
ALL_PPI_1 <- c(ALL_PPI_data[,1])
ALL_PPI_2 <- c(ALL_PPI_data[,2])
ALL_PPI_Proteins <- c(ALL_PPI_1,ALL_PPI_2)
Numbers_PPI_Proteins <- unique(ALL_PPI_Proteins)

#不同的共表达矩阵,两个跑一个就好了
#1、基于差异基因构建的共表达矩阵,pearson去相关性，相关性取0.5
co_expression_data <- read.table("F:/xuejie_data/Differential_Relation_rnaseq0.5-co_expression.txt",sep = " ")
#2.基于MIC前30%构建的共表达矩阵,MIC取相关性,相关性取0.5
co_expression_data <- read.table("F:/xuejie_data/MIC_Relation_rnaseq0.5-co_expression.txt",sep = " ")
#3、基于MIC前30%构建的共表达矩阵，pearson取相关性，相关性取0.5
co_expression_data <- read.table("F:/xuejie_data/Relation_rnaseq0.5-co_expression.txt",sep = " ")
#4、基于差异基因构建的共表达矩阵,MIC取相关性，相关性取0.5
co_expression_data <- read.table("F:/xuejie_data/0.5_Relation_rnaseq_Differexpression_MIC_Select.txt",sep = " ")

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

write.table(Intersection_Proteins,file = "E:/zhuliu_File/R Code_PPI/DATA/Mid_data/Intersection_Proteins.txt",quote = FALSE,sep = " ",col.names = FALSE,row.names = FALSE)
Intersection_Proteins <- read.csv("E:/zhuliu_File/R Code_PPI/DATA/Mid_data/Intersection_Proteins.txt",header = FALSE)


#对ALL_PPI和co_expression进行共同蛋白质筛选
PPI_final_1 <- as.data.frame(matrix(nrow=0,ncol=0))
PPI_final_2 <- as.data.frame(matrix(nrow=0,ncol=0))
PPI_final_1 <- subset(ALL_PPI_data,ALL_PPI_data[,1] %in% Intersection_Proteins[,1])
PPI_final_2 <- subset(PPI_final_1,PPI_final_1[,2] %in% Intersection_Proteins[,1])
write.table(PPI_final_2,file = "E:/zhuliu_File/R Code_PPI/DATA/PPI_2021_9_4.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)


co_final_1 <- as.data.frame(matrix(nrow=0,ncol=0))
co_final_2 <- as.data.frame(matrix(nrow=0,ncol=0))
co_final_1 <- subset(co_expression_data,co_expression_data[,1] %in% Intersection_Proteins[,1])
co_final_2 <- subset(co_final_1,co_final_1[,2] %in% Intersection_Proteins[,1])
write.table(co_final_2,file = "E:/zhuliu_File/R Code_PPI/DATA/expression_2021_9_4.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)




#判断有多少种子结点在这两个网络中
library(dplyr)
rm(list=ls())
setwd("E:/zhuliu_File/R Code_PPI")

All_PPI <- read.csv('DATA/PPI_2021_9_4.txt',sep = " ",header = FALSE)
co_expression_data <- read.csv('DATA/expression_2021_9_4.txt',sep = " ",header = FALSE)

Yi_Protein_data <- read.csv('DATA/Protein_Gene/nai_brought_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')

#All_PPI <- STRING_data[,1:2]
#All_PPI <- PPI_table

Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% All_PPI[,1] | Yi_Protein_data[,1] %in% All_PPI[,2])
Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% co_expression_data[,1] | Select_Seed_Proteins_1[,1] %in% co_expression_data[,2])

#Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% newdata_3[,1] | Yi_Protein_data[,1] %in% newdata_3[,2])
#Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% RicePPINET_data[,1] | Select_Seed_Proteins_1[,1] %in% RicePPINET_data[,2])
#Select_Seed_Proteins_3 <- subset(Select_Seed_Proteins_2,Select_Seed_Proteins_1[,1] %in% STRING_data[,1] | Select_Seed_Proteins_1[,1] %in% STRING_data[,2])

write.table(Select_Seed_Proteins_2,file = "DATA/Intersection_Seed_Proteins.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)






###############################################
#从表达矩阵出发，建立共表达网络
#1、从差异表达筛选后的表达矩阵出发，此数据分析0.5、0.6、0.7阈值下建立的共表达网络
rm(list = ls())
library(Hmisc)#加载包
#差异表达基因
mydata <- read.csv("F:/xuejie_data/Rna-seq_Differential_expression_Select.csv")
#MIC前30%，做pearson测试
#mydata <- read.csv("F:/xuejie_data/1-rnaseq/drought/myrice_drought_rnamic30.csv")
rownames(mydata) <- mydata$X
mydata <- mydata[,-(1)]
mydata <- t(mydata)
res <- cor(mydata)
round(res, 2)#保留两位小数

library(WGCNA)
probes <- colnames(res)
dimnames(res) <- list(probes,probes)
cyt <- exportNetworkToCytoscape(res,
                                edgeFile = paste("exprMAt",".edge.txt",sep = ""),
                                nodeFile = paste("exprMAt",".node.txt",sep = ""),
                                weighted = TRUE,threshold = 0,
                                nodeNames = probes)
woyaode <- cyt[["edgeData"]]
woyaode <- woyaode[,1:3]
woyaode <- subset(woyaode,woyaode$weight>0.6)
#write.table(woyaode,file = "F:/xuejie_data/MIC30_pearson_Relation_rnaseq0.6-co_expression.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)

write.table(woyaode,file = "F:/xuejie_data/Differential_Relation_rnaseq0.5-co_expression.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)



#2、从MIC前30%构建的表达矩阵出发，用MIC值来确定共表达网络
rm(list = ls())
library(minerva)
library(Hmisc)#加载包
mydata <- read.csv("F:/xuejie_data/1-rnaseq/drought/myrice_drought_rnamic30.csv")
rownames(mydata) <- mydata$X
mydata <- mydata[,-(1:2)]
MIC30_matrix <- mine(mydata)

MIC30_matrix_m <- as.matrix(MIC30_matrix)

MIC30_matrix_M <- matrix(0,nrow = 29521,ncol = 1)
MIC30_matrix_M <- MIC30_matrix_m[1,]#取出MIC
MIC30_matrix_M_MIC <- as.matrix(MIC30_matrix_M[["MIC"]])
rownames(MIC30_matrix_M_MIC) <- colnames(mydata)


write.csv(MIC30_matrix_M_MIC,file = "F:/xuejie_data/Relation_rnaseq_Differexpression_MIC_Select.csv",quote = FALSE,col.names = TRUE,row.names = TRUE)
write.csv(MIC30_matrix_M_MIC,file = "F:/xuejie_data/Relation_rnaseq_MIC_Select-co_expression.csv",quote = FALSE,col.names = TRUE,row.names = TRUE)
MIC30_matrix_M_MIC <- read.csv("F:/xuejie_data/Relation_rnaseq_MIC_Select-co_expression.csv")

MIC30_matrix_M_MIC <- read.csv("F:/xuejie_data/Relation_rnaseq_Differexpression_MIC_Select.csv")

library(WGCNA)
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
write.table(woyaode,file = "F:/xuejie_data/MIC_Relation_rnaseq0.4-co_expression.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)




write.table(woyaode,file = "F:/xuejie_data/0.5_Relation_rnaseq_Differexpression_MIC_Select.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)





#############################################################
#找四个不同网络中共同包含的种子节点，用来作为测试集，测试预测效果
rm(list = ls())
MIC_Relation_expression_data <- read.csv("F:/xuejie_data/MIC_Relation_rnaseq0.6-co_expression.txt",sep = " ")
MIC_Relation_expression_data_1 <- c(MIC_Relation_expression_data[,1])
MIC_Relation_expression_data_2 <- c(MIC_Relation_expression_data[,2])
MIC_Relation_expression_data_Proteins <- c(MIC_Relation_expression_data_1,MIC_Relation_expression_data_2)
Numbers_MIC_Relation_expression_data_Proteins <- unique(MIC_Relation_expression_data_Proteins)

Diff_expression_data <- read.table("F:/xuejie_data/Differential_Relation_rnaseq0.5-co_expression.txt",sep = " ")
Diff_expression_data_1 <- c(Diff_expression_data[,1])
Diff_expression_data_2 <- c(Diff_expression_data[,2])
Diff_expression_data_Proteins <- c(Diff_expression_data_1,Diff_expression_data_2)
Numbers_Diff_expression_data_Proteins <- unique(Diff_expression_data_Proteins)


ALL_PPI_data <- read.table("E:/zhuliu_File/R Code_PPI/DATA/ALL_PPI_Data.txt",sep = " ")
ALL_PPI_data_1 <- c(ALL_PPI_data[,1])
ALL_PPI_data_2 <- c(ALL_PPI_data[,2])
ALL_PPI_data_Proteins <- c(ALL_PPI_data_1,ALL_PPI_data_2)
Numbers_ALL_PPI_data_Proteins <- unique(ALL_PPI_data_Proteins)

STRING_data <- read.table("E:/zhuliu_File/R Code_PPI/DATA/STRING_PPI/STRING_0.txt",sep = " ")
STRING_data_1 <- c(STRING_data[,1])
STRING_data_2 <- c(STRING_data[,2])
STRING_data_Proteins <- c(STRING_data_1,STRING_data_2)
Numbers_STRING_data_Proteins <- unique(STRING_data_Proteins)


Yi_Protein_data <- read.csv('DATA/Protein_Gene/nai_brought_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')


Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,1] | Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,2])
Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,1] | Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,2])
Select_Seed_Proteins_3 <- subset(Select_Seed_Proteins_2,Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,1] | Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,2])
Select_Seed_Proteins_4 <- subset(Select_Seed_Proteins_3,Select_Seed_Proteins_3[,1] %in% STRING_data[,1] | Select_Seed_Proteins_3[,1] %in% STRING_data[,2])

write.table(Select_Seed_Proteins_4,file = "E:/zhuliu_File/R Code_PPI/DATA/2-Test_Seeds_brought_MIC0.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)




#############################################################
#找四个不同网络中共同包含的种子节点，用来作为测试集，测试预测效果
rm(list = ls())
MIC_Relation_expression_data <- read.csv("F:/xuejie_data/Relation_rnaseq0.5-co_expression.txt",sep = " ")
MIC_Relation_expression_data_1 <- c(MIC_Relation_expression_data[,1])
MIC_Relation_expression_data_2 <- c(MIC_Relation_expression_data[,2])
MIC_Relation_expression_data_Proteins <- c(MIC_Relation_expression_data_1,MIC_Relation_expression_data_2)
Numbers_MIC_Relation_expression_data_Proteins <- unique(MIC_Relation_expression_data_Proteins)

Diff_expression_data <- read.table("F:/xuejie_data/Differential_Relation_rnaseq0.5-co_expression.txt",sep = " ")
Diff_expression_data_1 <- c(Diff_expression_data[,1])
Diff_expression_data_2 <- c(Diff_expression_data[,2])
Diff_expression_data_Proteins <- c(Diff_expression_data_1,Diff_expression_data_2)
Numbers_Diff_expression_data_Proteins <- unique(Diff_expression_data_Proteins)


ALL_PPI_data <- read.table("E:/zhuliu_File/R Code_PPI/DATA/1-Threshold Proteins/ALL_PPI_Data_Thre_0.3.txt",sep = " ")
ALL_PPI_data_1 <- c(ALL_PPI_data[,1])
ALL_PPI_data_2 <- c(ALL_PPI_data[,2])
ALL_PPI_data_Proteins <- c(ALL_PPI_data_1,ALL_PPI_data_2)
Numbers_ALL_PPI_data_Proteins <- unique(ALL_PPI_data_Proteins)

STRING_data <- read.table("E:/zhuliu_File/R Code_PPI/DATA/1-Threshold Proteins/STRING0.3.txt",sep = " ")
STRING_data_1 <- c(STRING_data[,1])
STRING_data_2 <- c(STRING_data[,2])
STRING_data_Proteins <- c(STRING_data_1,STRING_data_2)
Numbers_STRING_data_Proteins <- unique(STRING_data_Proteins)


Yi_Protein_data <- read.csv('E:/zhuliu_File/R Code_PPI/DATA/Protein_Gene/nai_brought_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')


Select_Seed_Proteins_1 <- subset(Yi_Protein_data,Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,1] | Yi_Protein_data[,1] %in% MIC_Relation_expression_data[,2])
Select_Seed_Proteins_2 <- subset(Select_Seed_Proteins_1,Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,1] | Select_Seed_Proteins_1[,1] %in% Diff_expression_data[,2])
Select_Seed_Proteins_3 <- subset(Select_Seed_Proteins_2,Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,1] | Select_Seed_Proteins_2[,1] %in% ALL_PPI_data[,2])
Select_Seed_Proteins_4 <- subset(Select_Seed_Proteins_3,Select_Seed_Proteins_3[,1] %in% STRING_data[,1] | Select_Seed_Proteins_3[,1] %in% STRING_data[,2])

write.table(Select_Seed_Proteins_4,file = "E:/zhuliu_File/R Code_PPI/DATA/1-Test_Seeds_brought_P0.3.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)

