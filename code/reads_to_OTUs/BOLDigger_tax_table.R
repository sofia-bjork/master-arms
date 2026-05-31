#!/usr/bin/env Rscript

library(readxl)
library(stringr)
library(dplyr)

boldigger_dir <- "output/COI/BOLDigger"
if(!dir.exists(boldigger_dir)) dir.create(boldigger_dir)

lulu_dir <- "output/COI/MOTU"

# Boldigger file
# Read the API corrected sheetand remove the ">" in MOTU IDs 
bold <- as.data.frame(read_xlsx(path = file.path(lulu_dir, "boldigger3_data", "COI_cluster_reps_lulu_curated_identification_result.xlsx"),
                                sheet = 1))

colnames(bold)[1] <- "MOTU"

## Generate final taxonomy
phylum <- ifelse(bold$pct_identity > 85, bold$phylum, NA)
class <- ifelse(bold$pct_identity > 85, bold$class, NA)
order <- ifelse(bold$pct_identity > 85, bold$order, NA)
family <- ifelse(bold$pct_identity > 90, bold$family, NA)
genus <- ifelse(bold$pct_identity > 95, bold$genus, NA)
species <- ifelse(bold$pct_identity > 98, bold$species, NA)

tax_table <- cbind(MOTU = bold$MOTU, phylum, class, order, family, genus, species)

# write final tax table to file
write.table(tax_table, file = file.path(boldigger_dir, "boldigger_tax_table.txt"),
            sep = "\t", row.names = FALSE)

# Read LULU curated MOTU count table.
motu_tab <- read.table(file = file.path(lulu_dir, "lulu_motu_table_COI.txt"),
                        header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Sort the count table based on the order in the tax table and combine both tables
motu_tab <- motu_tab[order(match(motu_tab[, 1], tax_table[, 1])), ]

motu_tax_tab <- as.data.frame(cbind(tax_table, motu_tab[, -1]), stringsAsFactors = FALSE)

# Write this MOTU table to file
write.table(motu_tax_tab, file = file.path(boldigger_dir, "COI_motu_table_tax_counts_species_fullname.txt"), 
            sep = "\t", row.names = F)

## Aggregate MOTUs with same species assignment, adding read counts of same-species MOTUs to the most abundant same-species MOTU. ##

# Sort motu_tax_tab by read abundance 

motu_tax_tab <- motu_tax_tab[order(rowSums(motu_tax_tab[,9:ncol(motu_tax_tab)]), decreasing = TRUE),]

# Create a mapping file for subsequent downstream analysis to track which MOTUs are going to be merged.

motu_map <- aggregate(MOTU ~ species, data = motu_tax_tab, paste, collapse = ",")
write.table(motu_map, file = file.path(boldigger_dir, "motu_map_identical_species.txt"),
            sep = "\t", row.names = FALSE)

# Aggregate rows by species entries and sum read counts.
tax_sum <- aggregate(.~ motu_tax_tab$species, data = motu_tax_tab[,9:ncol(motu_tax_tab)], sum)

# Combine this table with the mapping table
motu_unique <- as.data.frame(cbind(motu_map, tax_sum[,-1]))

# Where same-species MOTUs have been merged, we only want to keep the MOTU ID of the most abundant one.
# Where ID entries with aggregated MOTUs exist, the first entry equals the most abundant MOTU so we delete the characters after (and incl.) the comma.

motu_unique$MOTU <- gsub(",.*", "", motu_unique$MOTU)

# Filter the previously generated motu_tax_tab table to only contain the MOTUs now listed in the motu_unique table.
motu_tax_tab_uniq <- motu_tax_tab[motu_tax_tab$MOTU %in% motu_unique$MOTU,]

# Order the motu_unique table to match the MOTU ID order of motu_tax_tab_uniq and create a final MOTU table with the taxonomy info of motu_tax_tab_uniq and the summed up read counts of motu_unique 
motu_unique <- motu_unique[order(match(motu_unique[,2], motu_tax_tab_uniq[,1])),]
final_motu_tab <- as.data.frame(cbind(motu_tax_tab_uniq[,1:8], motu_unique[,-c(1,2)]))

# Combine it with the entries that were not part of the aggregation procedure, e.g. MOTUs with no species level assignment
final_motu_tab <- rbind(final_motu_tab, motu_tax_tab[is.na(motu_tax_tab$species),])

# Write count table to file for phyloseq processing
write.table(final_motu_tab[, c(1,8:ncol(final_motu_tab))],
            file = file.path(boldigger_dir, "COI_motu_count_table_merged_species.txt"),
            sep = "\t", row.names = FALSE)

# Write taxonomy table to file for phyloseq processing
write.table(final_motu_tab[,1:7],
            file = file.path(boldigger_dir, "COI_motu_tax_table_merged_species.txt"),
            sep = "\t", row.names = FALSE)
