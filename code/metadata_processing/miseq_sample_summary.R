#!/usr/bin/env Rscript

# this script downloads the metadata files needed from the ARMS-MBON
# GitHub workspace and exports them to the right directory

# from the command line, one can specify what ObservatoryID or UnitID to look at
# also specify Preservative and Fraction of interest

library(argparse)
library(dplyr)

# specify meta and generated meta directories
meta_dir <- "metadata/meta"
gen_meta <- "metadata/generated_meta"

# specify (and create) meta directory
if(!dir.exists(meta_dir)) dir.create(meta_dir)

# specify (and create) output generated meta directory
if(!dir.exists(gen_meta)) dir.create(gen_meta)

# download metadata from github
combin_omics <- read.csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/refs/heads/main/qualitycontrolled_data/combined/combined_OmicsData.csv")
demult_omics <- read.csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/refs/heads/main/qualitycontrolled_data/combined/demultiplexing_details_OmicsData.csv")
sample_event <- read.csv("https://raw.githubusercontent.com/arms-mbon/data_workspace/refs/heads/main/qualitycontrolled_data/combined/combined_SamplingEventData.csv")

# save metadata files to the data/meta directory
write.csv(combin_omics, file = file.path(meta_dir, "combined_OmicsData.csv"),
          row.names = FALSE, col.names = TRUE, sep = ",")
write.csv(demult_omics, file = file.path(meta_dir, "demultiplexing_details_OmicsData.csv"),
          row.names = FALSE, col.names = TRUE, sep = ",")
write.csv(sample_event, file = file.path(meta_dir, "combined_SamplingEventData.csv"),
          row.names = FALSE, col.names = TRUE, sep = ",")

# Initialize command line argument parser
parser <- ArgumentParser(description = 'Accessing specified ENA sequences from metadata available at the ARMS-MBON GitHub.')

# All of your command line arguments
parser$add_argument('-o', '--observatory', metavar = 'ObservatoryID', type = 'character', required = FALSE, help = 'Specify ARMS ObservatoryID of interest.')
parser$add_argument('-p', '--preservative', metavar = 'Preservative', type = 'character', required = FALSE, help = 'Specify Preservative if needed.')
parser$add_argument('-u', '--unit', metavar = 'UnitID', type = 'character', required = FALSE, help = 'Specify UnitID of interest.')
parser$add_argument('-f', '--fraction', metavar = 'Fraction', type = 'character', required = FALSE, help = 'Specify Fraction of interest.')

# Parse the arguments
args <- parser$parse_args()

# Access the arguments
observ_id <- args$observatory
pres_spec <- args$preservative
unit_id <- args$unit
fraction <- args$fraction

# extract columns of interest from metadata files
combin_omics <- select(combin_omics, MaterialSampleID, Gene_COI,
                       Gene_COI_negative_control, Gene_COI_demultiplexed,
                       Gene_COI_comment, SequencingRunComment)

dim(combin_omics)

demult_omics <- select(demult_omics, MaterialSampleID, Sequencing_batch, Gene_COI)
dim(demult_omics)

sample_event <- select(sample_event, MaterialSampleID, ObservatoryID,
                       UnitID, Fraction, Filter, Preservative,
                       DateDeployed, DateCollected)
dim(sample_event)

# merge the three data bases based on MaterialSampleID
intermediate_merge <- merge(combin_omics, sample_event, by = "MaterialSampleID")
metadata_sum <- merge(intermediate_merge, demult_omics, by = "Gene_COI")

# if observ_id has an assigned character from the -o flag
if (!is.null(observ_id)) {
  metadata_sum <- metadata_sum[metadata_sum$ObservatoryID == observ_id, ]
}

# if pres_spec has an assigned character from the -p flag
if (!is.null(pres_spec)) {
  metadata_sum <- metadata_sum[metadata_sum$Preservative == pres_spec, ]
}

# if unit_id has an assigned character from the -u flag
if (!is.null(unit_id)) {
  metadata_sum <- metadata_sum[metadata_sum$UnitID == unit_id, ]
}

# if fraction has an assigned character from the -f flag
if (!is.null(fraction)) {
  metadata_sum <- metadata_sum[metadata_sum$Fraction == fraction, ]
}

# merge identical rows
metadata_sum <- unique(metadata_sum)

# export the summary file as csv
write.csv(metadata_sum, file = file.path(gen_meta, "COI_demultiplexed_summary.csv"),
          row.names = FALSE)
