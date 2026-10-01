## Script to generate a sample map between a vcf and the expression matrix generated from phaserPopPipe phaser_expr_matrix.py
## It's a little different for the dGTEx data, since the sample names are NOT the same between the bed file and the vcf

## Load libaries
library(data.table)
library(readxl)
library(dplyr)
library(tidyverse)

## Read in arguments 
args <- commandArgs(trailingOnly = TRUE)
samples <- args[1]
outpath <- args[2]
tissue <- unlist(strsplit(outpath, "/"))[2]

## Read in the samples FROM THE VCF
samples <- read.delim(samples, header=F)

## Now we need to bring back the metadata 
meta <- read_xlsx("/gpfs/commons/datasets/controlled/dGTEx/prerelease/dgtex_250728/dGTEx_Analysis_2025-07-21_v1/RNAseq_files/dGTEx_Analysis_2025-07-21_v1_RNA_Metadata_subset.xlsx") |> 
  as.data.frame()
meta$`Tissue Site` <- gsub(" ", "_", meta$`Tissue Site`)
## For the samples that we have in a given tissue, get their extended name from the bed file 
sub <- meta[meta$`Tissue Site`==tissue,]

map <- 
  sub |> 
  dplyr::select(`Donor ID`, `Sample ID (Data Release Official ID)`)
colnames(map) <- c("vcf_sample", "bed_sample")

## Make a sample map that has the IDs of the vcf and bed file sample
## Since we use the same names across these files, it's just a dataframe with two identical columns
# map <- data.frame(vcf_sample = samples$V1, 
#                   bed_sample = samples$V1)

## Write the output
write.table(map, outpath, row.names = FALSE, sep = "\t", quote = F)

