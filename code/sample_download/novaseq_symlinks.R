#!/usr/bin/env Rscript

library(dplyr)

# specify metadata dir for sample translation file
meta_dir <- "/cfs/klemming/home/b/bjorkso/Private/temporal-arms-main/metadata/generated_meta"

# specify project directory to symlink to 
fastq_dir <- "/cfs/klemming/home/b/bjorkso/Private/temporal-arms-main/output/novaseq_fastq"
if(!dir.exists(fastq_dir)) dir.create(fastq_dir)

# specify project directory to symlink from 
source_dir <- "/cfs/klemming/projects/supr/naiss2025-23-46/ARMS_Genoscope_data_full/projet_DBB/ARMS-Macro/"
pcr1_dir <- file.path(source_dir, "NEGATIVE_CONTROLS", "pcr1")
pcr1_files <- list.files(pcr1_dir)
pcr2_dir <- file.path(source_dir, "NEGATIVE_CONTROLS", "pcr2")
pcr2_files <- list.files(pcr2_dir)
code_dir <- file.path(source_dir, "COI", "COI_primers_m1COIintF_jgHCO2198")
code_files <- list.files(code_dir)


# read sample summary files from samples deployed in 2021 and 2022
sum_2122 <- read.csv("metadata/generated_meta/koster_coi_2021_2022.csv",
                         sep = ",", header = TRUE)
sum_2223 <- read.csv("metadata/generated_meta/koster_coi_2022_2023.csv",
                         sep = ",", header = TRUE)

# merge sample summary files
sample_sum <- bind_rows(sum_2122, sum_2223)

# gather all PCR negative control sample codes
pcr_control1 <- unique(sample_sum$PCR_negative_control_Code_1)
pcr_control2 <- unique(sample_sum$PCR_negative_control_Code_2)

# divide sample code into separate parts
code_split <- do.call(rbind, strsplit(sample_sum$Code, "_"))
code_split <- cbind(code_split, do.call(rbind, strsplit(code_split[, 4], "\\.")))

# recreate sample names (first two parts of Code strings)
sample_names <- paste(code_split[,1], code_split[,2], sep = "_")
sample_names <- unique(sample_names)
# extract sequencing batch and merge duplicates
sample_batch <- unique(code_split[, 5])
# assign a batch number to each sequencing run 
sample_batch <- cbind(c(1:length(sample_batch)), sample_batch)
colnames(sample_batch) <- c("Batch_X", "seq_run_code")

# function to direct target directory from filename
symlink_target_function <- function(file){
    for (i in 1:nrow(sample_batch)){
        if (grepl(sample_batch[i, 2], file) == TRUE){
            target_dir <- file.path(fastq_dir, paste0("Batch_", sample_batch[i, 1]))
            if(!dir.exists(target_dir)) dir.create(target_dir)
        }
    }

    if (is.null(target_dir)) {
        stop(message("No target dir found for file:", file))
    }

    return(target_dir)
}

# symlink pcr1 samples 
for (file in pcr_control1){
    pcr1_match <- grep(file, pcr1_files, value = TRUE)
    if (length(pcr1_match) > 1){
        for (file in pcr1_match){
            file.symlink(file.path(pcr1_dir, file), symlink_target_function(file))
        }
    }
}

# symlink pcr2 samples 
for (file in pcr_control2){
    pcr2_match <- grep(file, pcr2_files, value = TRUE)
    if (length(pcr2_match) > 1){
        for (file in pcr2_match){
            file.symlink(file.path(pcr2_dir, file), symlink_target_function(file))
        }
    }
}

# create list of code_files subdirectories to search code samples in
batch_list <- list()
for (i in 1:length(sample_batch[,2])){
    match <- grep(sample_batch[i, 2], code_files, value = TRUE)
    if (length(match) > 0) {
        batch_list <- c(batch_list, match)
    }
}

# search for sample names in each subdirectory in batch_list
for (seq_run in batch_list){
    subdir_path <- file.path(code_dir, seq_run)
    subdir_files <- list.files(subdir_path)
    for (name in sample_names){
            code_match <- grep(name, subdir_files, value = TRUE)
            if (length(code_match) > 1){
            for (file in code_match){
            file.symlink(file.path(subdir_path, file), symlink_target_function(file))
            }
            }
}
}

# export batch translations to metadata directory 
write.csv(sample_batch, file = file.path(meta_dir, "Sample_name_translations.csv"), 
          row.names = FALSE, quote = FALSE)

