################################################################################
#实现PPI三种不同数据库数据合并，阈值筛选
################################################################################
library(dplyr)
rm(list=ls())
setwd("E:/zhuliu_File/R Code_PPI")
PRIN_data <- read.csv('DATA/PRIN_PPI/PRIN_data_deal.csv')
RicePPINET_data <- read.csv('DATA/RicePPINET_PPI/RicePPINET_data_deal.csv')
STRING_data <- read.table('DATA/STRING_PPI/STRING_data_deal.txt',header = TRUE)


PRIN_data[,3] <- abs(PRIN_data[,3])

MAX_score <- max(STRING_data[,3])

STRING_data[,3] <- STRING_data[,3]/MAX_score

STRING_data <- data.frame(STRING_data)

colnames(PRIN_data) <- c('V1','V2','V3')
colnames(RicePPINET_data) <- c('V1','V2','V3')
colnames(STRING_data) <- c('V1','V2','V3')
################################################################################

Yi_Protein_data <- read.csv('DATA/Protein_Gene/nai_brought_genes.csv',header = TRUE)
colnames(Yi_Protein_data) <- c('V1','V2')
threshold_1 <- length(Yi_Protein_data[,1])
for (i in 0:10) {
  
  threshold <- i/10

  newdata_1 <- subset(PRIN_data,V3>threshold)
  newdata_2 <- subset(RicePPINET_data,V3>threshold)
  newdata_3 <- subset(STRING_data,V3>threshold)

  ALL_data_deal <- rbind(newdata_1[,1:2],newdata_2[,1:2],newdata_3[,1:2])

  All_PPI <- unique(ALL_data_deal)

  ################################################################################
  Yi_Protein_data_1 <- Yi_Protein_data$V1


  number_yizhi <- sum(Yi_Protein_data_1 %in% All_PPI[,1] | Yi_Protein_data_1 %in% All_PPI[,2])

  yizhi_protein_ratio <- number_yizhi/length(Yi_Protein_data[,1])

  
  #Calculate the ratio of seed nodes in the protein interaction network
  zhi <- c(i/10,yizhi_protein_ratio)
  zhi <- t(data.frame(zhi))
  hautu <- as.data.frame(matrix(nrow=0,ncol=0))
  if(i==0){
    hautu_1 <- rbind(hautu,zhi)
  }
  else{
    hautu_1 <- rbind(hautu_1,zhi)
  }
}

plot(hautu_1[,1],hautu_1[,2],type = 'b',xlab = 'threshold',ylab = 'percentage',main = 'Percentage of seed nodes under the threshold')

###########################################################################################################
threshold <- 0

newdata_1 <- subset(PRIN_data,V3>threshold)
newdata_2 <- subset(RicePPINET_data,V3>threshold)
newdata_3 <- subset(STRING_data,V3>threshold)


ALL_data_deal_1 <- rbind(PRIN_data[,1:2],RicePPINET_data[,1:2],STRING_data[,1:2])

#ALL_data_deal <- rbind(newdata_1[,1:2],newdata_2[,1:2],newdata_3[,1:2])

All_PPI <- unique(ALL_data_deal)

#write.table(newdata_3, file = 'DATA/STRING_0.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
#write.table(All_PPI, file = 'DATA/1-Threshold Proteins/ALL_PPI_Data_Thre_0.2.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
write.table(All_PPI, file = 'DATA/2-ALL Proteins/ALL_PPI_Data.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
  

