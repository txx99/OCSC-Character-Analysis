# OCSC-Character-Analysis
Tumor samples were collected from 4 patients with high-grade serous ovarian cancer (HGSOC). RNA-sequencing was performed on the day of biopsy collection (day 0) and again after 6 days of in vitro culturing.
Here, we investigate how well the in vitro samples retain the stem-like characteristics seen in the initial samples.
To do this, we perform {edgeR} DEA and {clusterProfiler} GSEA, with visualisations using a scatter plot, volcano plot, bar plot, and dot plot.

Raw Data: Read_Count_Data

Code: Preprocessing.R > blocked_DEA.R > GSEA.R

Outputs: Clean_Data > DEA_results.csv


Workflow: Preprocessing.R > Clean_Data/clean_all_data.csv | blocked_DEA.R > DEA_results.csv | GSEA.R
