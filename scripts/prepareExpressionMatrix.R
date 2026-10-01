## This script looks through the samplesheet for each tissue, grabs the phaser_gene_ae.txt files for those samples, and creates a directory that contains the necessary files for the expression matrix script
## The sample sheet path and the path to where the files are located should be in the config files

## Load libraries
library(data.table)
library(tidyverse)
library(dplyr)

## Set file paths and arguments 
args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
samplesheet <- args[2]
directory <- args[3]
ref <- args[4]
output <- args[5]

## Read in the sample sheet and subset it for tissue of interest 
samples <- fread(samplesheet) |> 
  filter(Tissue == tissue) |> 
  pull(Sample)

## For this tissue, get all of the gene_ae.txt files: 
files <- list.files(paste0(directory, tissue), full.names = TRUE)

## Select the samples of interest from this list of files: 
data <- files[unlist(lapply(samples, grep, files))] |> 
  lapply(fread)
samples <- files[unlist(lapply(samples, grep, files))] |> 
  basename()

## Get a list of all genes from all of the samples 
total <- do.call(rbind, data)
genes <- unique(total$name)

## Subset the reference for the genes in the above file: 
ref <- fread(ref) |> 
  filter(V4 %in% genes)

## down the line it'll cause issues if there's duplicate genes--remove genes that appear more than once in the reference 
remove <- names(which(table(ref$V4)>1))
ref <- ref[!ref$V4 %in% remove,]
## remove these genes from the gene list too--
genes <- unique(ref$V4)
## make sure the order is identical: 
identical(ref$V4, genes)

## write this output (THIS NEEDS TO BE IN THE SAME ORDER AS THE GENES )
write.table(ref, paste0(output, "input/", tissue, "/ref/", tissue, "_genes.bed"), row.names = FALSE, sep = "\t", quote = F, col.names = FALSE)


## For each file, leftjoin it with the reference 
for (i in 1:length(data)){
  print(i)
  merged <- left_join(ref, data[[i]], by = c("V4" = "name"), relationship = "many-to-many") |> 
    rename(name = V4) 
  merged <- merged[match(genes, merged$name),]
  merged <- merged[,c(1:4,17:24)] |> 
    mutate_if(is.numeric,coalesce,0)
  bam <- unique(merged$bam)
  bam <- bam[grep("final", bam)]
  merged$bam <- unlist(strsplit(samples[i], "_combined_ASE.txt"))
  merged[is.na(merged$variants),]$variants <- ""
  colnames(merged) <- c("contig", "start", "stop", "name", "aCount", "bCount", "totalCount", "log2_aFC", "n_variants", "variants", "gw_phased", "bam")
  merged <- merged[match(genes, merged$name),]
  sample <- unlist(strsplit(samples[i], "_combined_ASE.txt")) ## this will need to be changed to whatever the naming scheme for the other dataset has 
  out <- paste0(output, "input/", tissue, "/ae/", sample, "_sorted_ASE.txt")
  write.table(merged, out, row.names = FALSE, sep = "\t", quote = FALSE)
}
