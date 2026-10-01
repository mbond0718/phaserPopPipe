### Helper script for the phaser pop pipeline to generate QTL pairs objects for all tissues: 
 
## Script to make a pairs file for phaserPOP for GTEx eQTLs 
files <- list.files("/gpfs/commons/datasets/controlled/GTEx/portal/data/v8/GTEx_Analysis_v8_eQTL/", "egenes", full.names = TRUE)
tissues <- basename(files)
tissues <- unlist(strsplit(tissues, ".v8.egenes.txt.gz"))
outpath <- "/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPopPipe/ref/pairs/"

for (i in 1:length(files)){
  print(i)
  data <- fread(files[i])
  varsplit <- unlist(strsplit(data$variant_id, "_"))
  data$chrom <- varsplit[seq(1, length(varsplit)*5, by = 5)] |> na.omit()
  data$pos <- varsplit[seq(2, length(varsplit)*5, by = 5)] |> na.omit()
  data$ref <- varsplit[seq(3, length(varsplit)*5, by = 5)] |> na.omit()
  data$alt <- varsplit[seq(4, length(varsplit)*5, by = 5)] |> na.omit()
  ## get the hgnc gene name instead of ensembl id 
  ref <- fread("/gpfs/commons/groups/lappalainen_lab/mbond/work/asePipe_V3/ref/gencode.v25.annotation.bed")
  ref <- ref[,c(4,13)]
  colnames(ref) <- c("name", "gene_id")
  data <- left_join(data, ref, by = "gene_id")
  
  sub <- data[,c("name", "variant_id", "chrom", "pos", "ref", "alt")]
  colnames(sub) <- c("gene_id", "var_id", "var_contig", "var_pos", "var_ref", "var_alt")
  
  
  write.table(sub, paste0(outpath, tissues[i], "_pairs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}





tmp <- fread("/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPop/output/cerebellum/test.gw_phased.bed.gz")
## try to remove that gene that has EC10 in the name
tmp <- tmp[-grep("EC10", tmp$name),]
## Remove samples where the gene appears multiple times 
table <- table(tmp$name)
remove <- names(which(table>1))
tmp <- tmp[!tmp$name %in% remove,]


write.table(tmp, "/gpfs/commons/groups/lappalainen_lab/mbond/work/phaserPop/output/cerebellum/test.gw_phased_filt.bed", row.names = FALSE, sep = "\t", quote = F)

