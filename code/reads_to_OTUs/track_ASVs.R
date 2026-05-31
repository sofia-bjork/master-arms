#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)

####### track fate of ASVs thorughout the pipeline ###########
# specify path to output dir
post_dir <- "output/final_taxonomy"
if(!dir.exists(post_dir)) dir.create(post_dir)

## Create final mapping file with all ASVs that have been placed in one MOTU throughout the pipeline ##

# read LULU and swarm output files and the mapping file of merging MOTUs with same species assignment 

lulu <- read.table("output/COI/MOTU/motu_map_lulu_COI.txt", sep = "\t", header = T)
swarm <- read.table("output/COI/swarm/output.txt", sep = "\t", stringsAsFactors = F) 
spec_merge <- read.table("output/COI/BOLDigger/motu_map_identical_species.txt", sep = "\t", header = T)
mapping <- read.table("output/COI/BOLDigger/motu_asv_mapping.txt", sep = "\t", header = T)

# Separate swarm table into columns by space
# Remove read count strings from ASV names (lapply and function necessary to do this for every entry in the table, not just specific columns)

swarm <- separate_wider_delim(swarm,V1, delim = " ", names_sep = "", too_few = "align_start")
swarm[] <- lapply(swarm, function(y) gsub("_.*","", y))

# Order entries of first column in swarm table based on first column in lulu table
swarm <- swarm[order(match(swarm$V11, lulu$X)), ]

# Merge swarm and lulu tables
swarm_lulu_map<-cbind(lulu[,4],swarm)
colnames(swarm_lulu_map)[1]<-"MOTU"

# Separate MOTU strings in spec_merge into separate columns
spec_merge <- separate_wider_delim(spec_merge, MOTU, delim = ",", names_sep = "", too_few = "align_start")

# Create ID column in spec_merge (= MOTU other MOTUs have been merged onto) and transform to long format
spec_merge <- cbind(spec_merge[ ,c(2, 2:ncol(spec_merge))])
colnames(spec_merge)[1]<-"ID"

spec_merge_long <- melt(setDT(spec_merge), id.vars = "ID", variable.name = "string")
spec_merge_long <- spec_merge_long[ , -2]
colnames(spec_merge_long)[2] <- "MOTU"
spec_merge_long <- spec_merge_long[!is.na(spec_merge_long$MOTU), ]

## Merge swarm_lulu_map and spec_merge_long

swarm_lulu_spec_merge_map <- as.data.frame(merge(swarm_lulu_map, spec_merge_long, by = "MOTU", all = TRUE))
swarm_lulu_spec_merge_map <- swarm_lulu_spec_merge_map %>% relocate(ID)

# Where ID is NA, fill in with entry of MOTU column. Then, remove MOTU column.

swarm_lulu_spec_merge_map$ID <- ifelse(is.na(swarm_lulu_spec_merge_map$ID), swarm_lulu_spec_merge_map$MOTU, swarm_lulu_spec_merge_map$ID)
swarm_lulu_spec_merge_map <- swarm_lulu_spec_merge_map[ , -2]

# Write all columns except ID column into one string, MOTU names separated by comma

swarm_lulu_spec_merge_map$MOTU_string <- apply(swarm_lulu_spec_merge_map[ , 2:ncol(swarm_lulu_spec_merge_map)], 1,paste, collapse = ",") 

# Remove ,NA strings (occurred during the previous step when empty columns were pasted together)
# Keep only ID column and the MOTU_string column

swarm_lulu_spec_merge_map$MOTU_string <- gsub(",NA.*","", swarm_lulu_spec_merge_map$MOTU_string)
swarm_lulu_spec_merge_map <- swarm_lulu_spec_merge_map[ , c(1,ncol(swarm_lulu_spec_merge_map))]

# Aggregate rows based on ID column

swarm_lulu_spec_merge_map <- aggregate(.~ ID, data = swarm_lulu_spec_merge_map, paste, collapse = ",")

# Map ASVs to MOTUs 

# motu_asv_mapping.txt stems from initial phyloseq processing performed previously 
mapping <- mapping[order(match(mapping[ , 1], swarm_lulu_spec_merge_map[ , 1])), ]

swarm_lulu_spec_merge_map <- cbind(mapping[, 2], swarm_lulu_spec_merge_map)
colnames(swarm_lulu_spec_merge_map)[1:2] <- c("MOTU", "ASV_representative")
swarm_lulu_spec_merge_map <- separate_wider_delim(swarm_lulu_spec_merge_map, MOTU_string, delim = ",", names_sep = "", too_few = "align_start")

# summarize the number of ASVs per OTU and save to swarm_lulu_spec_merge_map
No_ASVs <- rowSums(!is.na(swarm_lulu_spec_merge_map[3:ncol(swarm_lulu_spec_merge_map)]))
swarm_lulu_spec_merge_map <- add_column(swarm_lulu_spec_merge_map, No_ASVs, .after = 2)

# save the merged ASV OTU file to csv
write.table(swarm_lulu_spec_merge_map, file = file.path(post_dir, "ASV_OTU_mapping.csv"),
            quote = FALSE, sep = ";", row.names = FALSE)

