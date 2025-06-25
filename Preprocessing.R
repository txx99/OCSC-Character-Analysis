# import required libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)

# instantiate read count tables 
A778_D0=read.table(file="1_CancerTranscriptomics/read_counts/A778_D0_gene.tsv", sep = '\t', header=TRUE)
A778_D6=read.table(file="1_CancerTranscriptomics/read_counts/A778_D6_gene.tsv", sep = '\t', header=TRUE)

A820_D0=read.table(file="1_CancerTranscriptomics/read_counts/A820_D0_gene.tsv", sep = '\t', header=TRUE)
A820_D6=read.table(file="1_CancerTranscriptomics/read_counts/A820_D6_gene.tsv", sep = '\t', header=TRUE)

A870_D0=read.table(file="1_CancerTranscriptomics/read_counts/A870_D0_gene.tsv", sep = '\t', header=TRUE)
A870_D0=read.table(file="1_CancerTranscriptomics/read_counts/A870_D6_gene.tsv", sep = '\t', header=TRUE)

A899_D0=read.table(file="1_CancerTranscriptomics/read_counts/A899_D0_gene.tsv", sep = '\t', header=TRUE)
A899_D0=read.table(file="1_CancerTranscriptomics/read_counts/A899_D6_gene.tsv", sep = '\t', header=TRUE)
