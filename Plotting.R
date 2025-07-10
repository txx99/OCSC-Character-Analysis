# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# -----library -------
# BiocManager::install("edgeR")
library(edgeR)
library(limma) #dependency of edgeR


# --- load cleaned data ----
patient_csv_files<-(list.files('./', "^clean_A", full.names = FALSE ))
all_data_file<- "clean_all_data.csv"

patient_dfs=list() # list of all clean patient dfs
for (i in patient_csv_files){
  len<- length(patient_dfs)
  name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) # name of df
  patient_dfs[[len+1]] <- assign(name, read.csv(i, row.names = 1)) #create independent 'name' df and append it to list
  patient_dfs[[len+1]] <- name # update element name in list 
  name<- get(name) 
  name[is.na(name)]<- 0
}

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
                          genes = "gene_name")
#------normalisation-----------------


#---------dispersion-------------------

# --------pairwise testing---------------