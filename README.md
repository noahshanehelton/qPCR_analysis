# qPCR_analysis

This tool will be used to analyze qPCR data in a rigorous and reproducible manner. 
Primary functionalities will include: 
**Primer Efficiency Calculations**, **Pfaffl Method Analysis**, and **Percent enrichment in polysome profiles**


Example data can be found in the data folder of this repository. Three datasets are included to represent want the different kinds of formats should be. The primer efficiency dataset is a csv file with the name "dilution_testdata.csv". The headers include Gene, Dilution, Ct1, Ct2, Ct3. The Pfaffl method dataset can be found in the "test_data.csv" file with headers: Gene, Condition, Ct1, Ct2, Ct3. The polysome profile dataset can be found in the file titled: polysome_profile_testdata.csv, with the headers including: "Gene	Fraction	Condition	Ct1	Ct2	Ct3" . A tutorial can be found in " Tutorial.ipynb " 
