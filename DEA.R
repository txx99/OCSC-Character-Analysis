# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR


# --- load cleaned data ----
patient_csv_files<-(list.files("^clean_A", full.names = FALSE ))
all_data_file<- "./Clean_Data/clean_all_data.csv"

patient_dfs=list() # list of all clean patient dfs
for (i in patient_csv_files){
  len<- length(patient_dfs)
  name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) # name of df
  patient_dfs[[len+1]] <- assign(name, read.csv(paste0('./Clean_Data/',i), row.names = 1)) #create independent 'name' df and append it to list
  patient_dfs[[len+1]] <- name # update element name in list 
  name<- get(name) 
  name[is.na(name)]<- 0 #assign NAs as 0s
}

# assign NAs as 0s
all_data<- read.csv(all_data_file, row.names = 1)
all_data[is.na(all_data)]<- 0

#duplicate gene names present so go by ensembl ids
# for (i in clean_patient_dfs){
#         i<- get(i)
#   #       print(any(duplicated(i$gene_name)))
#   #       print(any(duplicated(i$ensembl_gene_rec)))
#   #     }
# }


#-----boxplots + scatterplots? ------

# individual d0 v d6 boxplots
for (id in patient_dfs){
  df<- get(id)
  boxplot(df[,2:3], main = paste0("Boxplot of ", id, " Expression Values"), col="lightblue", outline=FALSE,
          ylab='Read Count')
}

boxplot(all_data[,2:9], main = paste0("Boxplot of All Expression Values"), col="lightblue", outline=FALSE,
          ylab='Read count')

# ----DGEList------
# create DGEList class object 
DGE_counts <- DGEList(counts = all_data,
                          genes = all_data$gene_name)


# add 'day' info as condition 
design_table<- as.data.frame(colnames(all_data[2:9]))
design_table<- rename(design_table, 'samples'='colnames(all_data[2:9])')
design_table$day<- as.factor(c(rep("day 0",4), rep("day 6",4)))

# Confirm order matches 
summary(colnames(all_data[2:9]) == design_table$samples)

# Add grouping information to DGEList object
DGE_counts$samples$group<- as.factor(design_table$day)

#------filter low expression genes---------
# filterByExpr keeps genes that have at least min.count ( = 10) reads in a worthwhile number samples
# will automatically look in your DGEList object for group info, etc
genes_to_keep <- filterByExpr(DGE_counts)
summary(genes_to_keep)

# filter data 
DGE_counts <- DGE_counts[genes_to_keep, keep.lib.sizes = FALSE] #recalculate library sizes
dim(DGE_counts)[1]==summary(genes_to_keep)[3]


relevant_genes<- list()
relevant_genes <- DGE_counts$genes %>% filter(gene_name == "ALDH1A1"| 
                                           gene_name =="CD44"| 
                                           # gene_name =="CD133"| 
                                           gene_name == "CD24"| 
                                           # gene_name =="CD117"| 
                                           gene_name == "EPCAM"| 
                                           gene_name ==  "THY1")
relevant_genes #all OCSC genes present 

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
DGE_counts<- estimateDisp(DGE_counts, design = model.matrix(~condition_))

# --------pairwise testing---------------
# Exact test measuring difference in mean of selected 'pair' of conditions
DGE_pairwise<- exactTest(DGE_counts, pair = c("day 0", "day 6")) #ensure these are saved as factors in DGEList object

# topTags() returns the top differentially expressed genes ('tags')
top_genes<- topTags(DGE_pairwise, n=20)


# ----find OCSC marker genes-------

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
# all OCSC have insignificant FDR
# CD44 top ranked of the five