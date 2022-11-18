# Rice_Genes_Prediction
Data of the Paper

The data of PPI download from three databases:      <br>
PRIN(http://bis.zju.edu.cn/prin/)： PRIN_data_deal.csv  <br>
RicePPINet(http://netbio.sjtu.edu.cn)： RicePPINET_data_deal.csv  <br>
STRING(https://cn.string-db.org/cgi/download?sessionId=bNjzN9zZTSdk)：   4530.protein.links.v11.5.txt.gz  <br>
Because the data is too large, you can download it from the website by yourself. The biological selection is Oryza sativa

Co_expression constructed by RNA-seq obtained from ncbi： tqnnew_drought.csv  <br>

The drought genes obtained from China Rice Data Center： nai_drought_genes.csv  <br>



The Code to deal with these data:  <br>
deal with the PPI data:   deal data.R    <br>
deal with intersection of two network:    Test_Code.R   <br>
LOOCV  and predict potential genes :     R packages_Test.R  <br>
