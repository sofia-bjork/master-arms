#!/usr/bin/env Rscript

library(dada2)

# this script will do chimera removal first on miseq samples, then on the novaseq samples
# the cleaned ASV tables are then merged to create one unified ASV data set 

# specify input directories
miseq_input <- "output/miseq_fastq"
novaseq_input <- "output/novaseq_fastq"

# specify output directories
miseq_dir <- file.path("output", "COI", "miseq_dada2")
novaseq_dir <- file.path("output", "COI", "novaseq_dada2")
merged_out <- "output/COI"

# for miseq: read COI sample summary and extract run numbers
coi_summary <- read.csv("metadata/generated_meta/COI_demultiplexed_summary.csv", 
                        sep = ",", header = TRUE)

# for novaseq: load sequencing batch numbers
batch_dir    <- file.path("output", "novaseq_fastq")
batch_list <- list.files(batch_dir, pattern = "Batch")




####################### miseq chimera removal #######################

# Load the sequence tables of the different sequence runs
coi_run <- unique(coi_summary$Sequencing_batch)

rds_list <- list()
for (run in coi_run) {
  run_file <- file.path(miseq_dir, paste("seqtab_Run", run, ".rds", sep = ""))
  rds_run <- readRDS(run_file)
  rds_list[[as.character(run)]] <- rds_run
}

# Merge sequence tables

merged <- do.call(mergeSequenceTables, unname(rds_list))
saveRDS(merged, file = file.path(miseq_dir, "merged_seqtab_coi.rds"))

# Keep sequence reads with a length of 310, 313 or 316 bp only.

seqtab.filtered <- merged[,nchar(colnames(merged)) %in% c(310, 313, 316)]
saveRDS(seqtab.filtered, file = file.path(miseq_dir, "seqtab_filtered_coi.rds"))

# Remove chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab.filtered, multithread=T, verbose=TRUE)

# Save sequence table with the non-chimeric sequences as RDS file:
saveRDS(seqtab.nochim, file = file.path(miseq_dir, "seqtab_nochim_coi.rds"))

# It is possible that a large fraction of the total number of UNIQUE SEQUENCES will be chimeras.
# However, this is usually not the case for the majority of the READS.
# Calculate percentage of the reads that were non-chimeric.

non_chimeric <- (sum(seqtab.nochim)/sum(merged))*100
message(non_chimeric, "% of reads were non-chimeric")

# Remove singletons from the non-chimeric ASVs

#Transform counts to numeric (as they will most likely be integers)
mode(seqtab.nochim) = "numeric"

# Subset columns with counts of > 1 and save to file
seqtab.nochim.nosingle <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
saveRDS(seqtab.nochim.nosingle, file = file.path(miseq_dir, "seqtab_nochim_nosingle_coi.rds"))

# Write a fasta file of the final, non-chimeric , non-singleton sequences with short >ASV... type headers

asv_seqs <- colnames(seqtab.nochim.nosingle)
asv_headers <- vector(dim(seqtab.nochim.nosingle)[2], mode = "character")
for (i in 1:dim(seqtab.nochim.nosingle)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(miseq_dir, "COI_nochim_nosingle_ASVs.fa"))

# Write an ASV count table of the final, non-chimeric, non-singleton sequences with short >ASV... type names

colnames(seqtab.nochim.nosingle) <- paste0("ASV", seq(ncol(seqtab.nochim.nosingle)))

ASV_counts<-t(seqtab.nochim.nosingle) # transposing table

write.table(ASV_counts, file = file.path(miseq_dir, "COI_ASV_counts_nosingle.txt"), sep = "\t", quote = F, col.names = NA)

## Track reads through the entire dada2 pipeline ##

# Track reads through the chimera and singleton removal step.

# Read non-chimeric and non-singleton table again, as it has been modified
seqtab.nochim.nosingle <- readRDS(file = file.path(miseq_dir, "seqtab_nochim_nosingle_coi.rds"))

track_nochim_nosingle <- cbind(rowSums(seqtab.filtered), rowSums(seqtab.nochim), rowSums(seqtab.nochim.nosingle))
colnames(track_nochim_nosingle) <- c("length_filt", "nonchim", "nosingle")

# Read tracking tables of the single runs and combine

track_list <- list()
for (run in coi_run) {
  track_file <- file.path(miseq_dir, paste("track_Run", run, ".txt", sep = ""))
  txt_track <- read.table(track_file, sep = "\t", header = T, row.names = 1)
  track_list[[as.character(run)]] <- txt_track
}

tracks <- do.call(rbind, unname(track_list))

# Combine all tracking tables

tracks <- tracks[order(match(rownames(tracks), rownames(track_nochim_nosingle))),]
track_coi <- cbind(tracks, track_nochim_nosingle)

# Calculate percentages for each step compared to input

track_coi$cutadapt_perc <- (track_coi$cutadapt / track_coi$input)*100
track_coi$filtered_perc <- (track_coi$filtered / track_coi$input)*100
track_coi$denoisedF_perc <- (track_coi$denoisedF / track_coi$input)*100
track_coi$denoisedR_perc <- (track_coi$denoisedR / track_coi$input)*100
track_coi$merged_perc <- (track_coi$merged / track_coi$input)*100
track_coi$length_filt_perc <- (track_coi$length / track_coi$input)*100
track_coi$nonchim_perc <- (track_coi$nonchim / track_coi$input)*100
track_coi$nosingle_perc <- (track_coi$nosingle / track_coi$input)*100

# Save final read tracking table to file

write.table(track_coi, file = file.path(merged_out, "track_miseq.txt"), sep = "\t", col.names = NA)





####################### novaseq chimera removal #######################

get_sample_name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
coi_run <- unname(sapply(batch_list, get_sample_name))

# Load the sequence tables of the different sequence runs

rds_list <- list()
for (run in coi_run) {
  run_file <- file.path(novaseq_dir, paste("seqtab_Batch", run, "_mod4.rds", sep = ""))
  rds_run <- readRDS(run_file)
  rds_list[[as.character(run)]] <- rds_run
}

# Merge sequence tables

if (length(rds_list) > 1) {
  # Merge sequence tables
  merged <- do.call(mergeSequenceTables, unname(rds_list))
  # Keep sequence reads with a length of 310, 313 or 316 bp only.
  seqtab.filtered <- merged[,nchar(colnames(merged)) %in% c(310, 313, 316)]
} else { 
  merged <- as.matrix(rds_list[[1]])
   # Keep sequence reads with a length of 310, 313 or 316 bp only.
  seqtab.filtered <- merged[,nchar(colnames(merged)) %in% c(310, 313, 316)]
}

saveRDS(merged, file = file.path(novaseq_dir, "merged_seqtab_coi.rds"))
saveRDS(seqtab.filtered, file = file.path(novaseq_dir, "seqtab_filtered_coi.rds"))

# Remove chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab.filtered, multithread=T, verbose=TRUE)

# Save sequence table with the non-chimeric sequences as RDS file:
saveRDS(seqtab.nochim, file = file.path(novaseq_dir, "seqtab_nochim_coi.rds"))

# It is possible that a large fraction of the total number of UNIQUE SEQUENCES will be chimeras.
# However, this is usually not the case for the majority of the READS.
# Calculate percentage of the reads that were non-chimeric.

non_chimeric <- (sum(seqtab.nochim)/sum(merged))*100
message(non_chimeric, "% of reads were non-chimeric")

# Remove singletons from the non-chimeric ASVs

#Transform counts to numeric (as they will most likely be integers)
mode(seqtab.nochim) = "numeric"

# Subset columns with counts of > 1 and save to file
seqtab.nochim.nosingle <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
saveRDS(seqtab.nochim.nosingle, file = file.path(novaseq_dir, "seqtab_nochim_nosingle_coi.rds"))

# Write a fasta file of the final, non-chimeric , non-singleton sequences with short >ASV... type headers

asv_seqs <- colnames(seqtab.nochim.nosingle)
asv_headers <- vector(dim(seqtab.nochim.nosingle)[2], mode = "character")
for (i in 1:dim(seqtab.nochim.nosingle)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(novaseq_dir, "COI_nochim_nosingle_ASVs.fa"))

# Write an ASV count table of the final, non-chimeric, non-singleton sequences with short >ASV... type names

colnames(seqtab.nochim.nosingle) <- paste0("ASV", seq(ncol(seqtab.nochim.nosingle)))

ASV_counts<-t(seqtab.nochim.nosingle) # transposing table

write.table(ASV_counts, file = file.path(novaseq_dir, "COI_ASV_counts_nosingle.txt"), sep = "\t", quote = F, col.names = NA)

## Track reads through the entire dada2 pipeline ##

# Track reads through the chimera and singleton removal step.

# Read non-chimeric and non-singleton table again, as it has been modified
seqtab.nochim.nosingle <- readRDS(file = file.path(novaseq_dir, "seqtab_nochim_nosingle_coi.rds"))

track_nochim_nosingle <- cbind(rowSums(seqtab.filtered), rowSums(seqtab.nochim), rowSums(seqtab.nochim.nosingle))
colnames(track_nochim_nosingle) <- c("length_filt", "nonchim", "nosingle")

# Read tracking tables of the single runs and combine

track_list <- list()
for (run in coi_run) {
  track_file <- file.path(novaseq_dir, paste("track_Batch", run, "_mod4.txt", sep = ""))
  txt_track <- read.table(track_file, sep = "\t", header = T, row.names = 1)
  track_list[[as.character(run)]] <- txt_track
}

tracks <- do.call(rbind, unname(track_list))

# Combine all tracking tables

tracks <- tracks[order(match(rownames(tracks), rownames(track_nochim_nosingle))),]
track_coi <- cbind(tracks, track_nochim_nosingle)

# Calculate percentages for each step compared to input

track_coi$cutadapt_perc <- (track_coi$cutadapt / track_coi$input)*100
track_coi$filtered_perc <- (track_coi$filtered / track_coi$input)*100
track_coi$denoisedF_perc <- (track_coi$denoisedF / track_coi$input)*100
track_coi$denoisedR_perc <- (track_coi$denoisedR / track_coi$input)*100
track_coi$merged_perc <- (track_coi$merged / track_coi$input)*100
track_coi$length_filt_perc <- (track_coi$length / track_coi$input)*100
track_coi$nonchim_perc <- (track_coi$nonchim / track_coi$input)*100
track_coi$nosingle_perc <- (track_coi$nosingle / track_coi$input)*100

# Save final read tracking table to file

write.table(track_coi, file = file.path(merged_out, "track_novaseq.txt"), sep = "\t", col.names = NA)

####################### merge ASV tables #######################

# read miseq and novaseq sequence tables

miseq_seqtab <- readRDS(file = file.path(miseq_dir, "seqtab_nochim_nosingle_coi.rds"))
novaseq_seqtab <- readRDS(file = file.path(novaseq_dir, "seqtab_nochim_nosingle_coi.rds"))

# merge sequence tables
merged_seqtab <- mergeSequenceTables(miseq_seqtab, novaseq_seqtab)

# Write a fasta file of the merged sequence tables with short >ASV... type headers
asv_seqs <- colnames(merged_seqtab)
asv_headers <- vector(dim(merged_seqtab)[2], mode = "character")
for (i in 1:dim(merged_seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(merged_out, "COI_nochim_nosingle_ASVs.fa"))

# Write an ASV count table of the merged sequence tables with short >ASV... type names
colnames(merged_seqtab) <- paste0("ASV", seq(ncol(merged_seqtab)))
ASV_counts<-t(merged_seqtab) # transposing table
write.table(ASV_counts, file = file.path(merged_out, "COI_ASV_counts_nosingle.txt"), sep = "\t", quote = F, col.names = NA)


