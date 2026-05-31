#!/usr/bin/env Rscript

library(dplyr)

# specify generated meta directory
gen_meta <- "metadata/generated_meta"

# load generated metadata file of samples of interest
# metadata_sum <- read.csv(file = file.path(gen_meta, "COI_demultiplexed_summary.csv"),
#                          sep = ",", header = TRUE)

metadata_sum <- read.csv(file = file.path(gen_meta, "COI_demultiplexed_summary.csv"), 
                         sep = ",", header = TRUE)

# extract Gene_COI and Sequencing_batch columns
metadata_sum <- select(metadata_sum, Gene_COI, Sequencing_batch, Gene_COI_negative_control)
# sort by Sequencing_batch from latest to first
metadata_sum <- metadata_sum[order(metadata_sum$Sequencing_batch,
                                   decreasing = TRUE), ]
# remove duplicate Gene_COI, keep the ones from the latest Sequencing_batch
metadata_sum <- metadata_sum[!duplicated(metadata_sum$Gene_COI),]
# sort by Sequencing_batch just in case
metadata_sum <- metadata_sum[order(metadata_sum$Sequencing_batch,
                                   decreasing = FALSE), ]

# divide the metadata_sum data frame into two fractions depending on
# number of characters in sample name
run_7digits <- metadata_sum[nchar(metadata_sum$Gene_COI) == 10, ]
run_8digits <- metadata_sum[nchar(metadata_sum$Gene_COI) == 11, ]

# extract the sequencing batches present in the chosen samples
# loop the sequencing numbers and pair to accession numbers
run_formatting <- function(parts) {
  run_numbers <- unique(parts$Sequencing_batch)

  # create an empty vector
  accessions <- c()

  for (n in run_numbers){
    accessions <- append(accessions, paste0(">Run_", n))

    samples <- parts$Gene_COI[parts$Sequencing_batch == as.integer(n)]

    control <- parts$Gene_COI_negative_control[parts$Sequencing_batch == as.integer(n)]
    control <- unique(control)

    accessions <- append(accessions, samples)
    accessions <- append(accessions, control)
    accessions <- append(accessions, "")
  }
  accessions <- accessions[-length(accessions)]
  return(accessions)
}

# run the function for both 7 digit and 8 digit accession numbers
file_7digits <- run_formatting(run_7digits)
file_8digits <- run_formatting(run_8digits)

write.table(file_7digits, file = file.path(gen_meta, "COI_ENA_accessions.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(file_8digits, file = file.path(gen_meta, "COI_ENA_accessions_8digits.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
