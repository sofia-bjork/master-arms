#!/usr/bin/env Rscript

library(devtools)
install_github("tobiasgf/lulu")

library(lulu)

# specify COI/MOTU directory

MOTU_dir <- "output/COI/MOTU"

# Read MOTU table

motu_table <- read.table(file = file.path(MOTU_dir, "motu_table_COI.txt"), sep = "\t", header = T, row.names = 1, as.is = T)

# Read match list

matchlist <- read.table(file = file.path(MOTU_dir, "match_list.txt"), sep="\t", header = F, as.is = T, stringsAsFactors = F)

# Run curation

curated_result <- lulu(motu_table, matchlist,minimum_match = 0.84,minimum_relative_cooccurence = 0.9)

# Write curated MOTU table to file

write.table(curated_result$curated_table, file = file.path(MOTU_dir, "lulu_motu_table_COI.txt"), sep = "\t", quote = F, col.names = NA)

# Write data on how MOTUs where mapped to file

write.table(curated_result$otu_map, file = file.path(MOTU_dir, "motu_map_lulu_COI.txt"), sep = "\t", quote = F, col.names = NA)

# Write headers of curated MOTUs to file to subset the corresponding representative sequences of the swarm cluster fasta file later on

lulu_curated_headers <- paste0(">",row.names(curated_result$curated_table))
write.table(lulu_curated_headers, file = file.path(MOTU_dir, "lulu_curated_headers.txt"), sep = "\t", row.names = F, quote = F, col.names = F)
