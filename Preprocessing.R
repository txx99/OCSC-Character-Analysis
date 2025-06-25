# import required libraries
library(dplyr)
library(tidyverse)
library(DESeq2)
library(pheatmap)

# instantiate read count tables 
A778_D0=read.table(file="Read_Count_Data/A778_D0_gene.tsv", sep = '\t', header=FALSE)
A778_D6=read.table(file="Read_Count_Data/A778_D6_gene.tsv", sep = '\t', header=FALSE)

A820_D0=read.table(file="Read_Count_Data/A820_D0_gene.tsv", sep = '\t', header=FALSE)
A820_D6=read.table(file="Read_Count_Data/A820_D6_gene.tsv", sep = '\t', header=FALSE)

A870_D0=read.table(file="Read_Count_Data/A870_D0_gene.tsv", sep = '\t', header=FALSE)
A870_D0=read.table(file="Read_Count_Data/A870_D6_gene.tsv", sep = '\t', header=FALSE)

A899_D0=read.table(file="Read_Count_Data/A899_D0_gene.tsv", sep = '\t', header=FALSE)
A899_D0=read.table(file="Read_Count_Data/A899_D6_gene.tsv", sep = '\t', header=FALSE)
