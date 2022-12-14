# Rice_Genes_Prediction
Data of the Paper

The data of PPI download from three databases:      <br>
PRIN(http://bis.zju.edu.cn/prin/)： PRIN_data_deal.csv  <br>
RicePPINet(http://netbio.sjtu.edu.cn)： RicePPINET_data_deal.csv  <br>
STRING(https://cn.string-db.org/cgi/download?sessionId=bNjzN9zZTSdk)：   4530.protein.links.v11.5.txt.gz  <br>
Because the data is too large, you can download it from the website by yourself. The biological selection is Oryza sativa

Co_expression constructed by RNA-seq obtained from ncbi： rice_drought_RNA_seq_data.csv  <br>

The drought genes obtained from China Rice Data Center： known_drought_stress_related_genes.csv  <br>

Differentially expressed genes： Rna-seq_Differential_expression_Select.csv    <br>

The top 30% MIC values genes: myrice_drought_rnamic30.csv    <br>

The Code to deal with these data:  <br>
deal with the PPI data:   Fusion_of_Multi-source_PPI_data.R    <br>
Code for constructing gene co-expression networks： Co_expression construction.R   <br>
deal with intersection of two network and buliding a test genes set:    Genes_Select.R   <br>
LOOCV and predict potential genes :     LOOCV_and_RWR.R  <br>
