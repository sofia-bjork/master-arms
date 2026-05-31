#!/usr/bin/env Rscript

### COI ###

output_dir <- "output/COI/blank_corr"
if(!dir.exists(output_dir)) dir.create(output_dir)

# read COI sample summary and extract negative control IDs for miseq samples
miseq_summary <- read.csv("metadata/generated_meta/COI_demultiplexed_summary.csv", 
                        sep = ",", header = TRUE)
miseq_negative <- unique(miseq_summary$Gene_COI_negative_control)

# read COI sample summary and extract negative control IDs for novaseq samples
novaseq_summary <- read.csv("metadata/generated_meta/koster_coi_2022_2023.csv", 
                        sep = ",", header = TRUE)
novaseq_negative <- c(unique(novaseq_summary$PCR_negative_control_Code_1), unique(novaseq_summary$PCR_negative_control_Code_2))

# merge COI sample summaries
coi_negative <- c(miseq_negative, novaseq_negative)

# Read ASV count table (output from dada2)

asv_table <- read.table("output/COI/COI_ASV_counts_nosingle.txt",
                        sep = "\t", header = T, row.names = 1,
                        as.is = T, check.names = F)

# Read list of ASVs which remained after NUMT removal

asv_list <- read.table("output/COI/pseudo/nonpseudo.combined.names.txt", sep =  "\t")

# Remove ">" in asv_list

asv_list[,1] <- gsub(">", "", asv_list[,1])

# Subset asv_table to non-numt ASVs

asv_table <- subset(asv_table, rownames(asv_table) %in% asv_list[,1])

# Subset rows where the ASV read count in the blank / negative samples  exceeds 10 % of an ASVs total read count 
# before combining the tables using rbind, add a column with ASV IDs. 
# This is necessary to stop R from introducing new rownames adding zeros to ASV names if duplicate ASVs exist in the newly created blank dataframes

blank_list <- list()

for (id in coi_negative) {
  subset_id <- asv_table[grep(id, colnames(asv_table))]

  negative_id <- subset(asv_table, subset_id > 0.1 * rowSums(asv_table))
  blank <- cbind(ASV_ID = rownames(negative_id), negative_id)

  blank_list[[id]] <- blank
}

# Combine blank dataframes
blanks <- do.call(rbind, unname(blank_list))

# Remove ASVs if there are duplicates in the "blanks" table (in case an ASV's read count exceeded 10 % of the total read count in more than one blank sample)

blanks <- blanks[!duplicated(blanks$ASV_ID), ]

# Remove the potential contaminant ASVs from the ASV table

asv_table <- cbind(ASV_ID = rownames(asv_table), asv_table)
asv_no_contams <- asv_table[!(asv_table$ASV_ID %in% blanks$ASV_ID),]

# filter the blank corrected ASVs 
# set reads representing less than 0.005% of the reads within a sample to 0
for (index in 2:ncol(asv_no_contams)){
  sample <- asv_no_contams[, index]
  tot_reads <- sum(sample)
  if (tot_reads == 0) next 

  fraction_vector <- sample/tot_reads
  asv_no_contams[, index][fraction_vector < 0.00005] <- 0
}

# identify samples with < 5 reads for an OTU and set those to 0
asv_no_contams[asv_no_contams < 5] <- 0

# remove OTUs with 0 reads
asv_no_contams <- asv_no_contams[rowSums(asv_no_contams[, -1]) > 0, ]

# Write table of potential contaminant ASVs and ASV count table devoid of contaminants

write.table(blanks, file = file.path(output_dir, "asv_contaminants_COI.txt"), sep="\t", row.names = F)

write.table(asv_no_contams, file = file.path(output_dir, "asv_no_contaminants_COI.txt"), sep="\t", row.names = F)

# Write headers of non-contaminant asvs to file to subset the corresponding sequences of non-NUMT fasta file later on

no_contam_headers <- paste0(">",asv_no_contams$ASV_ID)
write.table(no_contam_headers, file = file.path(output_dir, "no_contam_headers_COI.txt"), sep="\t", row.names = F,quote=F,col.names = F)
