# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

# --- load cleaned data ----

all_data_file<- "./Clean_Data/clean_all_data.csv"

# assign NAs as 0s
all_data<- read.csv(all_data_file)
all_data[is.na(all_data)]<- 0

#-----boxplots + scatterplots? ------
all_data_long <- all_data %>% 
  pivot_longer(!c(ensembl_gene_rec, gene_name), 
               names_to = c('patient', 'day'), 
               names_sep = '_',
               values_to = 'read_count')

all_data_wide <- all_data_long %>% 
  pivot_wider(names_from = day, 
              values_from = read_count)

all_data_wide %>%
  mutate(D0 = as.numeric(unlist(D0)),
         D6 = as.numeric(unlist(D6)),
         patient = as.factor(patient)) %>%
  ggplot(aes(x = D0, 
             y = D6, 
             colour = patient)) +
  geom_point(alpha = 0.5) +
  labs(x = 'D0 Read Count', 
       y = 'D6 Read Count') +
  theme_bw()

# boxplot(all_data[,2:9], main = paste0("Boxplot of All Expression Values"), col="lightblue", outline=FALSE,
#         ylab='Read count')

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

DEA_results <- topTags(lrt, n = Inf)$table
write.csv(DEA_results, 'DEA_results.csv', row.names = TRUE)

# select top 20 genes
top_genes<- topTags(lrt, n=20)


# ---- find OCSC marker genes-------

all_genes<- topTags(lrt, n=Inf)
OCSC_markers<- list()
OCSC_markers <- all_genes$table %>% filter(genes == "ALDH1A1"| 
                                              genes =="CD44"| 
                                              genes == "CD24"| 
                                              genes == "EPCAM"| 
                                              genes ==  "THY1")
OCSC_markers # stats for all OCSC genes
rm(all_genes)
# all OCSC have insignificant FDR
# CD44 top ranked of the five


OCSC_markers %>% ggplot(aes(x = genes, 
                            y = logFC)) + 
  geom_col() + 
  geom_text(aes(y = 2, 
                label = round(FDR, digits = 4))) +
  labs(x = 'gene', 
       y = 'logFC Value') +
  theme_classic()


OCSC_markers %>% EnhancedVolcano(lab = 'genes', 
                                 x = 'logFC' , 
                                 y = 'FDR',
                                 title = '',
                                 xlab = 'Log2 FC',
                                 ylab = 'FDR (p-adj)',
                                 pCutoff = 0.05, 
                                 FCcutoff = 2)


top_genes$table %>% EnhancedVolcano(lab = top_genes$table$genes, 
                                 x = 'logFC' , 
                                 y = 'FDR',
                                 title = '',
                                 xlab = 'Log2 FC',
                                 ylab = 'FDR (p-adj)',
                                 pCutoff = 0.05, 
                                 FCcutoff = 2)










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
