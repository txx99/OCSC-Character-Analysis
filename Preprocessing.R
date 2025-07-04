# BINF5503 Final Project
# Authors: Usman Ahmed, Tazmeen Gill, Noha El-Haj
# Dataset: 1 - Cancer Transcriptomics

# ----- import required libraries -----
library(tidyverse)

# ----- import data -----
# patient 1
A778_D0 <- read_tsv(file = "Read_Count_Data/A778_D0_gene.tsv", 
                    col_names = c("ensembl_gene_rec",
                                  "gene_name",
                                  "read_count"),
                    na = "0")
A778_D6 <- read_tsv(file = "Read_Count_Data/A778_D6_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")

# patient 2
A820_D0 <- read_tsv(file = "Read_Count_Data/A820_D0_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")
A820_D6 <- read_tsv(file = "Read_Count_Data/A820_D6_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")

# patient 3
A870_D0 <- read_tsv(file = "Read_Count_Data/A870_D0_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")
A870_D6 <- read_tsv(file = "Read_Count_Data/A870_D6_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")

# patient 4
A899_D0 <- read_tsv(file = "Read_Count_Data/A899_D0_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")
A899_D6 <- read_tsv(file = "Read_Count_Data/A899_D6_gene.tsv", 
                    col_names = c("ensembl_gene_rec", 
                                  "gene_name", 
                                  "read_count"), 
                    na = "0")

# ----- create individual df for each patient -----
A778_df <- merge(A778_D0, 
                 A778_D6, 
                 by = c("ensembl_gene_rec", "gene_name"))
colnames(A778_df) <- c("ensembl_gene_rec", "gene_name", "d0_rc", "d6_rc") # rc = read count

A820_df <- merge(A820_D0, 
                 A820_D6, 
                 by = c("ensembl_gene_rec", "gene_name"))
colnames(A820_df) <- c("ensembl_gene_rec", "gene_name", "d0_rc", "d6_rc")

A870_df <- merge(A870_D0, 
                 A870_D6, 
                 by = c("ensembl_gene_rec", "gene_name"))
colnames(A870_df) <- c("ensembl_gene_rec", "gene_name", "d0_rc", "d6_rc")

A899_df <- merge(A899_D0, 
                 A899_D6, 
                 by = c("ensembl_gene_rec", "gene_name"))
colnames(A899_df) <- c("ensembl_gene_rec", "gene_name", "d0_rc", "d6_rc")

# create list of all patients for ease of function application
patients_list <- mget(ls(pattern = "^A.*_df$"))

# ----- data exploration -----
lapply(patients_list, head)
lapply(patients_list, summary) 
lapply(patients_list, str)
lapply(patients_list[[1]], typeof)
lapply(patients_list[[1]], class)

# check for OCSC genes in dataset
relevant_genes <- list()
for (i in 1:length(patients_list)) {
  relevant_genes[[i]] <- patients_list[[i]] %>% filter(gene_name == "ALDH1A1"| 
                                                         gene_name =="CD44"| 
                                                         # gene_name =="CD133"| 
                                                         gene_name == "CD24"| 
                                                         # gene_name =="CD117"| 
                                                         gene_name == "EPCAM"| 
                                                         gene_name ==  "THY1")
}






