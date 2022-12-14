################################################################################
#Function is used for the fusion of multi-source protein interaction data.
################################################################################

library(dplyr)
rm(list=ls())

#The data includes three columns, Genes A, Genes B, and a score between A and B.
PRIN_data <- read.csv('PRIN_data_deal.csv')
RicePPINET_data <- read.csv('RicePPINET_data_deal.csv')
STRING_data <- read.table('STRING_data_deal.txt',header = TRUE)


PRIN_data[,3] <- abs(PRIN_data[,3])
MAX_score <- max(STRING_data[,3])
STRING_data[,3] <- STRING_data[,3]/MAX_score
STRING_data <- data.frame(STRING_data)

colnames(PRIN_data) <- c('V1','V2','V3')
colnames(RicePPINET_data) <- c('V1','V2','V3')
colnames(STRING_data) <- c('V1','V2','V3')

#the relation between A and B is considered, if the score more than 0.3. 
threshold <- 0.3

newdata_1 <- subset(PRIN_data,V3>threshold)
newdata_2 <- subset(RicePPINET_data,V3>threshold)
newdata_3 <- subset(STRING_data,V3>threshold)


ALL_data_deal_1 <- rbind(PRIN_data[,1:2],RicePPINET_data[,1:2],STRING_data[,1:2])

All_PPI <- unique(ALL_data_deal)

write.table(All_PPI, file = 'ALL_PPI_Data_0.3.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
  

