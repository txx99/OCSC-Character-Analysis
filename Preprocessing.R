# BINF5503 Final Project
# Authors: Usman Ahmed, Tazmeen Gill, Noha El-Haj
# Dataset: 1 - Cancer Transcriptomics

# ----- import required libraries -----
library(dplyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)

# ----- import data -----
A778_D0 <- read_tsv(file = "Read_Count_Data/A778_D0_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))
A778_D6 <- read_tsv(file = "Read_Count_Data/A778_D6_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))

A820_D0 <- read_tsv(file = "Read_Count_Data/A820_D0_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))
A820_D6 <- read_tsv(file = "Read_Count_Data/A820_D6_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))

A870_D0 <- read_tsv(file = "Read_Count_Data/A870_D0_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))
A870_D6 <- read_tsv(file = "Read_Count_Data/A870_D6_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))

A899_D0 <- read_tsv(file = "Read_Count_Data/A899_D0_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))
A899_D6 <- read_tsv(file = "Read_Count_Data/A899_D6_gene.tsv", 
                    col_names = c("Ensembl Gene Record", 
                                  "Common Gene Name", 
                                  "Read Count"))

# ----- merge data into one dataframe -----

