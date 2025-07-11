# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR


# --- load cleaned data ----

setwd('./Clean_Data')
getwd()
patient_csv_files<-(list.files(pattern="^clean_A", full.names = FALSE ))
setwd('..')

patient_dfs=list() # list of all clean patient dfs
for (i in patient_csv_files){
  len<- length(patient_dfs)
  name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) # name of df
  patient_dfs[[len+1]] <- assign(name, read.csv(paste0('./Clean_Data/',i), row.names = 1)) #create independent 'name' df and append it to list
  patient_dfs[[len+1]] <- name # update element name in list 
  name_df<- get(name)
  name_df<- rename(name_df, "day 0"="d0_rc", 'day 6'='d6_rc')
  name_df[is.na(name_df)]<- 0
  assign(name, name_df) #assign NAs as 0s
}


#-----boxplots + scatterplots? ------

# individual d0 v d6 boxplots
for (id in patient_dfs){
  df<- get(id)
  boxplot(df[,2:3], main = paste0("Boxplot of ", id, " Expression Values"), col="lightblue", outline=FALSE,
          ylab='Read Count')
}

# ----DGEList------

# create DGEList class object 
for (i in patient_dfs){
  design_table<- as.data.frame(colnames(df[2:3]))
  design_table<- rename(design_table, 'samples'='colnames(df[2:3])')
  design_table$day<- as.factor(c("day 0","day 6"))
  
  df<- get(i)
  DGE_counts <- DGEList(counts = df,
                      genes = df$gene_name)
  DGE_counts$samples$group<- as.factor(design_table$day)
  summary(colnames(df[2:3]) == design_table$samples)
  DGE_counts$samples$group<- as.factor(design_table$day)
  
  genes_to_keep <- filterByExpr(DGE_counts)
  summary(genes_to_keep)
  DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
  
  condition_<- design_table$day
  DGE_counts<- estimateDisp(DGE_counts, design = model.matrix(~condition_))
  DGE_pairwise<- exactTest(DGE_counts, pair = c("day 0", "day 6")) #ensure these are saved as factors in DGEList object
  
  top_genes<- topTags(DGE_pairwise, n=20)
  most_genes<- topTags(DGE_pairwise, n=15000)
  OCSC_markers<- list()
  OCSC_markers <- most_genes$table %>% filter(gene_name == "ALDH1A1"| 
                                                gene_name =="CD44"| 
                                                # gene_name =="CD133"| 
                                                gene_name == "CD24"| 
                                                # gene_name =="CD117"| 
                                                gene_name == "EPCAM"| 
                                                gene_name ==  "THY1")
  OCSC_markers # stats for all OCSC genes
  
  assign(paste0(i, '_DGE_pairwise'), DGE_pairwise)
  assign(paste0(i, '_top_genes'), top_genes)
  assign(paste0(i, '_OCSC_markers'), OCSC_markers)
  
}


# 
# # Confirm order matches 
# summary(colnames(df[2:3]) == design_table$samples)
# 
# # Add grouping information to DGEList object
# DGE_counts$samples$group<- as.factor(design_table$day)
# 
# #------filter low expression genes---------
# # filterByExpr keeps genes that have at least min.count ( = 10) reads in a worthwhile number samples
# # will automatically look in your DGEList object for group info, etc
# genes_to_keep <- filterByExpr(DGE_counts)
# summary(genes_to_keep)
# 
# # filter data 
# DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
# dim(DGE_counts)[1]==summary(genes_to_keep)[3]
# 
# 
# relevant_genes<- list()
# relevant_genes <- DGE_counts$genes %>% filter(gene_name == "ALDH1A1"| 
#                                                 gene_name =="CD44"| 
#                                                 # gene_name =="CD133"| 
#                                                 gene_name == "CD24"| 
#                                                 # gene_name =="CD117"| 
#                                                 gene_name == "EPCAM"| 
#                                                 gene_name ==  "THY1")
# relevant_genes #all OCSC genes present 
# 
# # ----------normalisation [SKIP?]-----------------
# # comparing gene expr bw samples, but different samples have different library sizes (total read count) --> can skew analysis
# # normalise so that each sample has relatively similar impact on DEA analysis, reduce bias for high expr genes
# DGE_counts$samples$norm.factors #currently all samples weighed equally
# 
# # calcNormFactors() normalises normalizes across all samples based on lib sizes (assumption: genes *not* differentially expressed)
# # by scaling to minimise LogFC bw samples of the same treatment
# DGE_counts <- calcNormFactors(DGE_counts)
# DGE_counts$samples$norm.factors # see normalisation factors have been adjusted
# 
# #---------dispersion-------------------
# # estimate gene dispersion == estimate relative variability of true expr bw replicates
# # create design matrix (good practice)
# condition_<- design_table$day
# DGE_counts<- estimateDisp(DGE_counts, design = model.matrix(~condition_))
# 
# # --------pairwise testing---------------
# # Exact test measuring difference in mean of selected 'pair' of conditions
# DGE_pairwise<- exactTest(DGE_counts, pair = c("day 0", "day 6")) #ensure these are saved as factors in DGEList object
# 
# # topTags() returns the top differentially expressed genes ('tags')
# top_genes<- topTags(DGE_pairwise, n=20)
# 
# 
# # ----find OCSC marker genes-------
# 
# most_genes<- topTags(DGE_pairwise, n=15000)
# OCSC_markers<- list()
# OCSC_markers <- most_genes$table %>% filter(gene_name == "ALDH1A1"| 
#                                               gene_name =="CD44"| 
#                                               # gene_name =="CD133"| 
#                                               gene_name == "CD24"| 
#                                               # gene_name =="CD117"| 
#                                               gene_name == "EPCAM"| 
#                                               gene_name ==  "THY1")
# OCSC_markers # stats for all OCSC genes
# # all OCSC have insignificant FDR
# # CD44 top ranked of the five