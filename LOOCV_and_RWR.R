##########################################################
#Author: Liu Zhu
#LOOCV_and_RWR: Function to 
#This section was used to do LOOCV testing,
#and RWR-M-based genetic prediction of drought stress in rice.
##########################################################


rm(list = ls())
#gc()
library(RandomWalkRestartMH)
library(base)
library(Matrix)
library(igraph)
library(Rcpp)
library(plyr)
library(ggplot2)

################################################################################
#Running RWR-M on three networks to evaluate the predicted performance.
#1.Constructing a PPI network on STRING database.
PPI_table <- read.table('STRING0.5.txt',header = TRUE)
PPI_table <- read.table('STRING_data_deal.txt',header = TRUE)
PPI_table <- PPI_table[,1:2]
PPI_table <- unique(PPI_table)
PPI_Network <- graph.data.frame(PPI_table,directed=FALSE)
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
PPI_Network <- create.multiplex(PPI_Network)

#2.Constructing a PPI network on multi-source data.
PPI_table <- read.csv('DALL_PPI_Data_Thre_0.3.txt',sep = " ")
PPI_table <- read.csv('ALL_PPI_Data.txt',sep = " ")
PPI_table <- PPI_table[,1:2]
PPI_table <- unique(PPI_table)
PPI_Network <- graph.data.frame(PPI_table,directed=FALSE)
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
PPI_Network <- create.multiplex(PPI_Network)

#3.Constructing a multiplex network on the data shared by PPI network and co-expression network.
PPI_table <- read.csv("PPI_2021_9_4.txt",header = FALSE,sep = " ")
co_expression_table <- read.csv("expression_2021_9_4.txt",header = FALSE,sep = " ")
co_expression_table <- co_expression_table[,1:2]

PPI_Network <- graph.data.frame(PPI_table,directed=FALSE)
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)

co_expression_Network <- graph.data.frame(co_expression_table,directed = FALSE)
co_expression_Network <- igraph::simplify(co_expression_Network,remove.multiple = TRUE, remove.loops = TRUE)

PPI_Network <- create.multiplex(PPI_Network,co_expression_Network)

################################################################################


Adjacenecy.matrix <- compute.adjacency.matrix(PPI_Network)
Adjacency.matrix.normalize <- normalize.multiplex.adjacency(Adjacenecy.matrix)


#Read in seed node set or test seed node set.
Seed_Nodes <- read.csv("Intersection_Seed_Proteins.txt",sep = " ",header = FALSE)
Seed.Nodes <- c(unique(Seed_Nodes[,1]))

Top_Random_walk_Results_Test<- as.data.frame(matrix(nrow=0,ncol=0))

#Two loop testing algorithms
for (item in 1:length(Seed.Nodes)) {
  Seed.Nodes.Deal <- Seed.Nodes[-item]
  Random_algorithm <- Random.Walk.Restart.Multiplex(Adjacency.matrix.normalize,PPI_Network,Seed.Nodes.Deal,r = 0.7)
  Random_walk_Results_Deal <- Random_algorithm[["RWRM_Results"]]
  
  One_Seed_Rank <- which(Seed.Nodes[item] == Random_walk_Results_Deal[,1])
  One_Seed_Rank_percentage <- One_Seed_Rank/(length(Random_walk_Results_Deal[,1])+length(Seed.Nodes))
  One_Seed_Rank_1 <- c(Seed.Nodes[item],One_Seed_Rank)
  One_Seed_Rank_percentage_1 <- c(Seed.Nodes[item],One_Seed_Rank_percentage)
  DATA_One_Rank <- as.data.frame(matrix(nrow=0,ncol=0))
  
  if(item == 1){
    LOOCV_Result <- rbind(DATA_One_Rank,One_Seed_Rank_1)
    LOOCV_Result_percentage <- rbind(DATA_One_Rank,One_Seed_Rank_percentage_1)
  }else{
    LOOCV_Result <- rbind(LOOCV_Result,One_Seed_Rank_1)
    LOOCV_Result_percentage <- rbind(LOOCV_Result_percentage,One_Seed_Rank_percentage_1)
  }
}
write.table(LOOCV_Result,file = 'Drought_related_LOOCV_Result.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
LOOCV_Result <- read.csv('Drought_related_LOOCV_Result.txt',sep = " ",header = FALSE)
colnames(LOOCV_Result) <- c("variable","value")
LOOCV_Result[,1] <- c("PPI")
LOOCV_Result_1 <- LOOCV_Result
LOOCV_Result_1 <- ddply(LOOCV_Result_1, .(variable), transform, ecd=ecdf(value)(value))


cdf <- ggplot(LOOCV_Result_1, aes(x=value)) + stat_ecdf(aes(colour=variable,linetype=variable),size = 1.75)+xlab("Rank") + 
              ylab("Cumulative Distribution") + theme_bw() +   theme(axis.text=element_text(size=16),
              axis.title=element_text(size=9,face="bold"), legend.title=element_blank(),
              legend.text = element_text(size = 12, face = "bold"), 
              legend.position=c(0.15,0.8), legend.key.size = unit(1.2, "cm"), 
              legend.background = element_rect(colour='black',size=0.5),
              panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
              panel.border=element_rect(size=0.75,colour='black')) 

cdf

cdf + coord_cartesian(xlim = c(1, 300), ylim = c(0,0.6))


#Runing RWR algorithm
Random_algorithm <- Random.Walk.Restart.Multiplex(Adjacency.matrix.normalize,PPI_Network,Seed.Nodes,r = 0.7)
Random_walk_Results <- Random_algorithm[["RWRM_Results"]]

#Retaining top n proteins according to LOOCV results
Random_walk_Results_Top_300 <- Random_walk_Results[1:300,]
write.table(Random_walk_Results_Top_300,file = 'DEG_M0.5_RWR_top300_Result.txt',sep = " ",row.names =FALSE, col.names =FALSE, quote =FALSE)
