# edgeR DEA reference:
# https://github.com/rnnh/bioinfo-notebook/blob/master/docs/DE_analysis_edgeR_script.md

# load cleaned data
patient_csv_files<-(list.files('./', "^clean_A", full.names = FALSE ))
day_csv_files<-(list.files('./', "^clean_d", full.names = FALSE ))

patient_dfs=list() # list of all clean patient dfs
for (i in patient_csv_files){
  len<- length(patient_dfs)
  name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) # name of df
  patient_dfs[[len+1]] <- assign(name, read.csv(i, row.names = 1)) #create independent 'name' df and append it to list
  patient_dfs[[len+1]] <- name # update element name in list 
}

day_dfs<- list()
for (i in day_csv_files){
  len<- length(day_dfs)
  name<- str_replace_all(i, c(".csv$"= "", "^clean_" = "")) 
  assign(name, read.csv(i, row.names = 1)) 
  day_dfs[[len+1]] <- name # update element name in list 
}

#duplicate gene names present so must go by ensembl ids
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

# all d0 vs all d6 boxplots
for (day in day_dfs){
  df<- get(day)
  boxplot(df[,2:5], main = paste0("Boxplot of ", day, " Expression Values"), col="lightblue", outline=FALSE,
          ylab='Read count')
}


#------normalisation-----------------


#---------dispersion-------------------

# --------pairwise testing---------------