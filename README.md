# Rice_Genes_Prediction
Data of the Paper

The data of PPI download from three databases: 
PRIN(http://bis.zju.edu.cn/prin/) : PRIN_data_deal.csv
RicePPINet(http://netbio.sjtu.edu.cn): RicePPINET_data_deal.csv
STRING(https://cn.string-db.org): STRING_data_deal.txt

Co_expression constructed by RNA-seq obtained from ncbi
MIC between Genes : MIC_Relation_rnaseq0.5-co_expression.txt
Pearson between Genes : Relation_rnaseq0.5-co_expression.txt

The drought genes obtained from China Rice Data Center : nai_drought_genes.csv



The Code to deal with these data:
deal data.R     #deal with the PPI data
Test_Code.R     #deal with intersection of two network
R packages_Test.R  #LOOCV  and predict potential genes
