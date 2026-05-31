#!/usr/bin/env Rscript

# For clustering ASVs into MOTUs using swarm, an ASV fasta with headers containing the total abundance of each is required.
# Here, we generate these new headers for the fasta files which will be used as input for swarm.

### COI ###

# Read the non-contaminant COI ASV count table 

ASV_counts <- read.table(file = "output/COI/blank_corr/asv_no_contaminants_COI.txt", sep = "\t", row.names = 1, header = T)

# Count total read abundances per ASV 

ASV_sums <-rowSums(ASV_counts)

# Create headers containing read counts

seqnames <- paste0(">", paste(rownames(ASV_counts), ASV_sums, sep="_"))

write.table(seqnames, "output/COI/ASV_dereplicated.txt", sep="\t", col.names = F, row.names = F, quote = F)

swarm_dir <- "output/COI/swarm"
if(!dir.exists(swarm_dir)) dir.create(swarm_dir)
