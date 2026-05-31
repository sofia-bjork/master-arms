#!/usr/bin/env Rscript

library(argparse)

# Initialize command line argument parser
parser <- ArgumentParser(description = 'COI FILTERING AND TRIMMING')

# All of your command line arguments
parser$add_argument('-d', '--directory', metavar = 'directory', type = 'character', required = TRUE, help = 'Specify cutadapt directory on machine.')

# Parse the arguments
args <- parser$parse_args()

# Access the arguments
cutadapt <- args$directory

# read step 2 in dx.doi.org/10.17504/protocols.io.n92ldmmmnl5b/v1 for more information about 
# which runs go through which function and why

# extract run information from summary csv
coi_summary <- read.csv("metadata/generated_meta/COI_demultiplexed_summary.csv",
                        sep = ",", header = TRUE)
coi_run <- unique(coi_summary$Sequencing_batch)

# for runs 1, 2, 5 and 6: 
# primer trimming and ASV inference with cutadapt and dada2 in R until read merging

source("code/reads_to_OTUs/miseq_filter_trim.R")
if (1 %in% coi_run == TRUE){
  seq_batch <- 1
  COI_trimming_filtering_func(seq_batch, cutadapt)
}

if (2 %in% coi_run == TRUE){
  seq_batch <- 2
  COI_trimming_filtering_func(seq_batch, cutadapt)
}

if (5 %in% coi_run == TRUE){
  seq_batch <- 5
  COI_trimming_filtering_func(seq_batch, cutadapt)
}

if (6 %in% coi_run == TRUE){
  seq_batch <- 6
  COI_trimming_filtering_func(seq_batch, cutadapt)
}

# for runs 3, 4 and 7:
# length filtering and ASV inference with cutadapt and dada2 until read merging
source("code/reads_to_OTUs/miseq_filter.R")
if (3 %in% coi_run == TRUE){
  seq_batch <- 3
  COI_filtering_func(seq_batch, cutadapt)
}

if (4 %in% coi_run == TRUE){
  seq_batch <- 4
  COI_filtering_func(seq_batch, cutadapt)
}

if (7 %in% coi_run == TRUE){
  seq_batch <- 7
  COI_filtering_func(seq_batch, cutadapt)
}
