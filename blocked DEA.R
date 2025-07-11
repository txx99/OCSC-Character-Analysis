# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR
library(tidyverse)


# --- load cleaned data ----

all_data_file<- "./Clean_Data/clean_all_data.csv"

# assign NAs as 0s
all_data<- read.csv(all_data_file, row.names = 1)
all_data[is.na(all_data)]<- 0

#-----boxplots + scatterplots? ------

boxplot(all_data[,2:9], main = paste0("Boxplot of All Expression Values"), col="lightblue", outline=FALSE,
        ylab='Read count')

# ----DGEList------
# create DGEList class object 
DGE_counts <- DGEList(counts = all_data,
                      genes = all_data$gene_name)
DGE_counts$genes$gene_name<-NULL

# add 'day' info as condition 
design_table<- data.frame(samples = colnames(all_data[2:9]))
design_table$patient<- as.factor(str_replace(design_table$samples, "_D[06]$",''))
design_table$day<- as.factor(c(rep("day 0",4), rep("day 6",4)))

# Confirm order matches 
summary(colnames(all_data[2:9]) == design_table$samples)

# Add grouping information to DGEList object
DGE_counts$samples$group<- as.factor(design_table$day)
DGE_counts$samples$group2<- as.factor(design_table$patient)


#------filter low expression genes---------
# filterByExpr keeps genes that have at least min.count ( = 10) reads in a worthwhile number samples
# will automatically look in your DGEList object for group info, etc
genes_to_keep <- filterByExpr(DGE_counts)
summary(genes_to_keep)

# filter data 
DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
dim(DGE_counts)[1]==summary(genes_to_keep)[3]

# ----------normalisation-----------------
# comparing gene expr bw samples, but different samples have different library sizes (total read count) --> can skew analysis
# normalise so that each sample has relatively similar impact on DEA analysis, reduce bias for high expr genes
DGE_counts$samples$norm.factors #currently all samples weighed equally

# calcNormFactors() normalises normalizes across all samples based on lib sizes (assumption: genes *not* differentially expressed)
# by scaling to minimise LogFC bw samples of the same treatment
DGE_counts <- calcNormFactors(DGE_counts)
DGE_counts$samples$norm.factors # see normalisation factors have been adjusted

#---------dispersion-------------------
# estimate gene dispersion == estimate relative variability of true expr bw replicates
# create design matrix (good practice)
condition_<- design_table$day
patient_<- design_table$patient

design_2 <- model.matrix(~patient_+condition_)
DGE_counts<- estimateDisp(DGE_counts, design = design_2)

# ---- glm test (blocking intrapatient) -------------
fit <- glmFit(DGE_counts, design_2)

# Likelihood ratio test
lrt <- glmLRT(fit, coef = "condition_day 6")
# design = matrix of what samples are being compared
# looking at: change in gene expression Day 0 to Day 6 after accounting for patient-specific differences

# select top 20 genes
top_genes<- topTags(lrt, n=20)


# ---- find OCSC marker genes-------

most_genes<- topTags(lrt, n=15000)
OCSC_markers<- list()
OCSC_markers <- most_genes$table %>% filter(genes == "ALDH1A1"| 
                                              genes =="CD44"| 
                                              genes == "CD24"| 
                                              genes == "EPCAM"| 
                                              genes ==  "THY1")
OCSC_markers # stats for all OCSC genes
rm(most_genes)
# all OCSC have insignificant FDR
# CD44 top ranked of the five






















# 
# 
# 
# 
# 
# 
# 
# setwd('./Clean_Data')
# getwd()
# patient_csv_files<-(list.files(pattern="^clean_A", full.names = FALSE ))
# setwd('..')
# 
# patient_dfs=list() # list of all clean patient dfs
# for (i in patient_csv_files){
#   len<- length(patient_dfs)
#   name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) # name of df
#   patient_dfs[[len+1]] <- assign(name, read.csv(paste0('./Clean_Data/',i), row.names = 1)) #create independent 'name' df and append it to list
#   patient_dfs[[len+1]] <- name # update element name in list 
#   name_df<- get(name)
#   name_df<- rename(name_df, "day 0"="d0_rc", 'day 6'='d6_rc')
#   name_df[is.na(name_df)]<- 0
#   assign(name, name_df) #assign NAs as 0s
# }
# 
# 
# #-----boxplots + scatterplots? ------
# 
# # individual d0 v d6 boxplots
# for (id in patient_dfs){
#   df<- get(id)
#   boxplot(df[,2:3], main = paste0("Boxplot of ", id, " Expression Values"), col="lightblue", outline=FALSE,
#           ylab='Read Count')
# }
# 
# # -----DGEList------
# # create DGEList class object
# 
# for (i in patient_dfs){
#   df<- get(i)
#   design_table<- as.data.frame(colnames(df[2:3]))
#   design_table<- rename(design_table, 'samples'='colnames(df[2:3])')
#   design_table$day<- as.factor(c("day 0","day 6"))
# 
# 
#   DGE_counts <- DGEList(counts = df,
#                       genes = df$gene_name)
#   DGE_counts$samples$group<- as.factor(design_table$day)
#   summary(colnames(df[2:3]) == design_table$samples)
#   DGE_counts$samples$group<- as.factor(design_table$day)
# 
#   genes_to_keep <- filterByExpr(DGE_counts)
#   summary(genes_to_keep)
#   DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
# 
#   condition_<- design_table$day
#   DGE_counts<- estimateDisp(DGE_counts, design = model.matrix(~condition_))
#   DGE_pairwise<- exactTest(DGE_counts, pair = c("day 0", "day 6")) #ensure these are saved as factors in DGEList object
# 
#   top_genes<- topTags(DGE_pairwise, n=20)
#   most_genes<- topTags(DGE_pairwise, n=15000)
#   OCSC_markers<- list()
#   OCSC_markers <- most_genes$table %>% filter(gene_name == "ALDH1A1"|
#                                                 gene_name =="CD44"|
#                                                 # gene_name =="CD133"|
#                                                 gene_name == "CD24"|
#                                                 # gene_name =="CD117"|
#                                                 gene_name == "EPCAM"|
#                                                 gene_name ==  "THY1")
#   OCSC_markers # stats for all OCSC genes
# 
#   assign(paste0(i, '_DGE_pairwise'), DGE_pairwise)
#   assign(paste0(i, '_top_genes'), top_genes)
#   assign(paste0(i, '_OCSC_markers'), OCSC_markers)
# 
# }
# 
