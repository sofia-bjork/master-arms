#!/usr/bin/env Rscript

library(phyloseq)
library(ggplot2)
library(data.table)


############# load relevant directories ##################

boldigger_dir <- "output/COI/BOLDigger"

# specify path to figures output
figures_dir <- "results/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)

# specify path to post-processing dir for worms taxon match file
post_dir <- "output/final_taxonomy"
if(!dir.exists(post_dir)) dir.create(post_dir)

# specify path to generated_meta directory (input and output directory)
meta_dir <- "metadata/generated_meta"


############# read necessary files for miseq summary file creation ##################

# read COI_demultiplexed_summary.csv for file names and sequencing run (miseq)
miseq_samples <- read.csv(file = file.path(meta_dir, "COI_demultiplexed_summary.csv"), sep = ",", header = TRUE)


############# read necessary files for novaseq summary file creation ##################

# read sample summary files from samples deployed in 2021 and 2022 and merge (novaseq)
sum_2122 <- read.csv(file = file.path(meta_dir, "koster_coi_2021_2022.csv"),
                         sep = ",", header = TRUE)
sum_2223 <- read.csv(file = file.path(meta_dir, "koster_coi_2022_2023.csv"),
                         sep = ",", header = TRUE)

# read sample name translations from symlink script (novaseq)
batch_translation <- read.csv(file = file.path(meta_dir, "Sample_name_translations.csv"),
                              sep = ",", header = TRUE)

# read track_COI.txt for sample names (novaseq)
track_COI <- read.csv("output/COI/track_novaseq.txt", sep = "\t", header = TRUE)
# read extracted_sample_names.csv for file names and sequencing run
extracted_samples <- read.csv("output/extracted_sample_names.csv", sep = ",", header = TRUE)



######### create COI sample summary files (novaseq) ###########

novaseq_sample_names <- track_COI$X
novaseq_sample_names <- cbind(sort(novaseq_sample_names), extracted_samples[order(extracted_samples$File_name), ])

# mapping matrix to transform dates to Month-Year format
month_map_novaseq <- data.frame(
  seq_run_code = c("HKYCCDRX3", "H3KVNDRX5", "HJ2LFDRX5", "H3WT3DRX7", "HLLHNDRX3", 
                   "H3L7NDRX5", "HJ5N5DRX5", "HCJF5DRX7", "HW5CMDRX3", "H55TTDRX5", "H3WLLDRX7"),
  Month_Year = c("Sep-23", "Apr-24", "Sep-24", "Oct-25", "Feb-24", 
                 "Apr-24", "Oct-24", "Nov-25", "Mar-24", "Jun-24", "Oct-25"),
  Seq_Date = c("230928", "240409", "240924", "251031", "240220",
               "240409", "241031", "251119", "240313", "240619", "251020"))

month_map_novaseq <- month_map_novaseq[order(month_map_novaseq$Seq_Date), ]

# create data frame for sample name (Sample), MaterialSampleID and month sequenced (Sequenced)
output <- data.frame(
  Sample = c(),
  MaterialSampleID = c(),
  Sequenced = c(),
  SeqMethod = c()
)
code <- append(sum_2122$Code, sum_2223$Code)
materialsampleid <- append(sum_2122$MaterialSample.ID..EMOBON., sum_2223$MaterialSample.ID..EMOBON.)
materialsampleid <- gsub("EMOBON", "ARMS", materialsampleid)
code_id_sums <- data.frame(code, materialsampleid)

for (sample in novaseq_sample_names[,1]){
  # extract row nr where the sample name matches the sample code
  row <- grep(sample, code_id_sums$code)
  # extract the MaterialSampleID from that same row
  id <- code_id_sums$materialsampleid[row]
  if (length(row) == 0){id <- "blank"}
  if (length(row) == 0){
    # look for sample name in novaseq_sample_names instead
    blank_row <- grep(sample, novaseq_sample_names$File_name)
    # extract sample codes to match the seq_run_code in month_map_novaseq
    code_split <- strsplit(novaseq_sample_names$File_name[blank_row], "_")
    code_split <- strsplit(code_split[[1]][5], "\\.")
    # extract the sequencing run code that matches the sample code
    month_row <- grep(code_split[[1]][1], month_map_novaseq$seq_run_code)
    month <- month_map_novaseq$Month_Year[month_row[1]]
  } else {
      # extract sample codes to match the seq_run_code in month_map_novaseq
    code_split <- strsplit(code_id_sums$code[row], "_")
    code_split <- strsplit(code_split[[1]][4], "\\.")
    # extract the sequencing run code that matches the sample code
    month_row <- grep(code_split[[1]][1], month_map_novaseq$seq_run_code)
    month <- month_map_novaseq$Month_Year[month_row[1]]
    }
  method <- "novaseq"
  # add all parts into the output data frame
  summary <- data.frame(
    Sample = sample, 
    MaterialSampleID = id,
    Sequenced = month,
    SeqMethod = method)
  output <- rbind(output, summary)
}

# strsplit the MaterialSampleIds to create Fraction, FractionGroup, 
# DeployedDate, RetrievedDate and RetrievedYear categories
output2 <- data.frame(
  Fraction = c(),
  FractionGroup = c(),
  DeployedDate = c(),
  RetrievedDate = c(),
  RetrievedYearYear = c(),
  ARMSUnit = c()
)
for (msID in output$MaterialSampleID){
  if (msID == "blank"){
    fraction_group <- "blank"
    fraction <- "blank"
    deployed_date <- "blank"
    retrieved_date <- "blank"
    retrieved_year <- "blank"
    arms_unit <- "blank"
  } else {
  msID_parts <- strsplit(msID, "_")
  fraction_group <- msID_parts[[1]][7]
  fraction <- gsub("[[:digit:]]", "", fraction_group)
  deployed_date <- msID_parts[[1]][5]
  retrieved_date <- msID_parts[[1]][6]
  retrieved_short <- substr(retrieved_date, 1, 2)
  retrieved_year <- paste0("20", retrieved_short)
  arms_unit  <- msID_parts[[1]][3]
  }
  summary2 <- data.frame(
  Fraction = fraction,
  FractionGroup = ifelse(fraction_group == "SF", "SF40", fraction_group),
  DeployedDate = deployed_date,
  RetrievedDate = retrieved_date,
  RetrievedYear = retrieved_year,
  ARMSUnit = arms_unit
  )
  output2 <- rbind(output2, summary2)
}

# final COI summary file
MOTUsampleCOI <- cbind(output, output2)
# save to generated_meta dir
# Write taxonomy table to file for phyloseq processing
write.table(MOTUsampleCOI,
            file = file.path(meta_dir, "novaseq_COI_sample_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)



######### create COI sample summary files (miseq) ###########

# create first miseq_COI_sample_data data frame
df1 <- data.frame(
  Sample = miseq_samples$Gene_COI,
  MaterialSampleID = miseq_samples$MaterialSampleID.x,
  Sequenced = miseq_samples$Sequencing_batch,
  SeqMethod = rep("miseq", nrow(miseq_samples))
)

# mapping matrix to transform sequencing batch to Month-Year format
month_map_miseq <- data.frame(
  seq_run = c(1, 2, 3, 4, 5, 6, 7),
  Month_Year = c("Jul-19", "Jan-20", "Sep-20", "Apr-21", "May-21", "Jan-22", "Aug-23"))

# replace sequencing run with Mon-YY month sequenced
for (run in month_map_miseq$seq_run){
  df1$Sequenced <- gsub(run, month_map_miseq$Month_Year[run], df1$Sequenced)
} 

df2 <- data.frame(
  Fraction <- c(),
  FractionGroup <- c(),
  DeployedDate <- c(),
  RetrievedDate <- c(),
  RetrievedYear <- c(),
  ARMSUnit <- c()
)
for (msID in df1$MaterialSampleID){
  # separate the parts of MaterialSampleID and make into column values
  msID_parts <- strsplit(msID, "_")
  fraction_group <- msID_parts[[1]][6]
  fraction <- gsub("[[:digit:]]", "", fraction_group)

  deployed_date_long <- msID_parts[[1]][4]
  deployed_date <- substr(deployed_date_long, 3, 8)

  retrieved_date_long <- msID_parts[[1]][5]
  retrieved_date <- substr(retrieved_date_long, 3, 8)

  retrieved_year <- substr(retrieved_date_long, 1, 4)

  arms_unit  <- msID_parts[[1]][3]
  # add all parts into the output data frame
  summary <- data.frame(
  Fraction = fraction,
  FractionGroup = fraction_group,
  DeployedDate = deployed_date,
  RetrievedDate = retrieved_date,
  RetrievedYear = retrieved_year,
  ARMSUnit = arms_unit
  )
  df2 <- rbind(df2, summary)
}

# merge the two metadata data frames
df3 <- cbind(df1, df2)

# change RetrievedYear for sample ERR4018451 from 2018 to 2019
df3$RetrievedYear[df3$Sample == "ERR4018451"] <- "2019"

# subset pcr negative control (blank) sample names and seq runs
sample_seqrun_subset <- data.frame(
  Sample <- miseq_samples$Gene_COI_negative_control,
  seq_run <- miseq_samples$Sequencing_batch
)
blank_samples <- unique(sample_seqrun_subset)

# create corresponding data frame for blanks
blanks <- data.frame(
  Sample = blank_samples$Sample,
  MaterialSampleID = rep("blank", nrow(blank_samples)),
  Sequenced = blank_samples$seq_run,
  SeqMethod = rep("miseq", nrow(blank_samples)),
  Fraction = rep("blank", nrow(blank_samples)),
  FractionGroup = rep("blank", nrow(blank_samples)),
  DeployedDate = rep("blank", nrow(blank_samples)),
  RetrievedDate = rep("blank", nrow(blank_samples)),
  RetrievedYear = rep("blank", nrow(blank_samples)),
  ARMSUnit = rep("blank", nrow(blank_samples))
) 

# change sequencing runs to sequencing date as done before
for (run in month_map_miseq$seq_run){
  blanks$Sequenced <- gsub(run, month_map_miseq$Month_Year[run], blanks$Sequenced)
} 

# merge blanks data frame with the rest of the samples and we have our metadata'
# final COI summary file
MOTUsampleCOI <- rbind(blanks, df3)
# save to generated meta directory
# Write taxonomy table to file for phyloseq processing
write.table(MOTUsampleCOI,
            file = file.path(meta_dir, "miseq_COI_sample_data.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)



######### create and clean phyloseq objects ###########

# load and merge  the two COI sample summary files

# Read sample metadata
miseq_sample_sum <- read.table(file = file.path(meta_dir, "miseq_COI_sample_data.txt"), 
                            header = T, row.names = 1, check.names = F, sep = "\t", strip.white = T)

novaseq_sample_sum <- read.table(file = file.path(meta_dir, "novaseq_COI_sample_data.txt"), 
                            header = T, row.names = 1, check.names = F, sep = "\t", strip.white = T)
# merge sample metadata
MOTUsampleCOI <- rbind(miseq_sample_sum, novaseq_sample_sum)

# Read MOTU counts 
MOTUcountsCOI <- read.table(file = file.path(boldigger_dir, "COI_motu_count_table_merged_species.txt"),
                            header = T, check.names = F, sep = "\t")

# Read MOTU taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)
MOTUtaxaCOI <- as.matrix(read.table(file = file.path(boldigger_dir, "COI_motu_tax_table_merged_species.txt"),
                                    header = T, check.names = F, sep = "\t"))

# Sort count table based on order in tax table (precautionary measure)
MOTUcountsCOI <- MOTUcountsCOI[order(match(MOTUcountsCOI[, 1], MOTUtaxaCOI[, 1])), ]

# MOTUs are still named "ASVxy". Replace them with MOTUxy and set them as rownames.

# First, create and write mapping file (MOTUxy = ASVyz). 
mapping <- cbind(MOTUcountsCOI[, 1], paste0("MOTU", seq(1:nrow(MOTUcountsCOI)))) 
colnames(mapping)[1:2] <- c("ASV","MOTU")
write.table(mapping, file = file.path(boldigger_dir, "motu_asv_mapping.txt"), sep = "\t", row.names = F)

rownames(MOTUcountsCOI) <- mapping[, 2]
MOTUcountsCOI[, 1] <- NULL

rownames(MOTUtaxaCOI) <- mapping[, 2]
MOTUtaxaCOI <- MOTUtaxaCOI[, -1] # setting it as NULL does not work for matrix object


# Create phyloseq object

psCOI <- phyloseq(otu_table(MOTUcountsCOI,taxa_are_rows = TRUE), sample_data(MOTUsampleCOI), tax_table(MOTUtaxaCOI))

## Get a quick overview of sequencing depth per sequencing run ##

# Make data.table for plot

read_sums <- data.table(as(sample_data(psCOI), "data.frame"),
                        TotalReads = sample_sums(psCOI), keep.rownames = TRUE)
setnames(read_sums, "rn", "SampleID")

# Violin plot based on sequencing events
sorted_months <- union(month_map_miseq$Month_Year, month_map_novaseq$Month_Year)

reads_plot <- ggplot(read_sums, aes(y=TotalReads,x=Sequenced,color=Sequenced)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2),size=1) +
  scale_x_discrete(limits=sorted_months)+
  theme_minimal() + 
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "#f5f3ed"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank()) +
  ggtitle("Sequencing Depth")

reads_plot

ggsave("COI_sequencing_depth_runs.png", plot = reads_plot, path = figures_dir,
       dpi = 400, width = 12, height = 12)

##

# Remove the negative controls
 
psCOI_cleaned <- subset_samples(psCOI,Fraction!="blank")


# make violin based on sequencing depths without blanks
read_sums <- data.table(as(sample_data(psCOI_cleaned), "data.frame"),
                        TotalReads = sample_sums(psCOI_cleaned), keep.rownames = TRUE)
setnames(read_sums, "rn", "SampleID")

reads_plot_no_blanks <- ggplot(read_sums, aes(y=TotalReads,x=Sequenced,color=Sequenced)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2),size=1) +
  scale_x_discrete(limits=sorted_months)+
theme_minimal() + 
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "#f5f3ed"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank()) +
  ggtitle("Sequencing Depth")

reads_plot_no_blanks

ggsave("COI_sequencing_depth_no_blanks.png", plot = reads_plot_no_blanks, path = figures_dir,
       dpi = 400, width = 12, height = 12)


# Remove the sediment samples and plankton samples (some sediment and plankton samples were sequenced as a trial during the initial phase of the ARMS program)

psCOI_cleaned <- subset_samples(psCOI_cleaned,Fraction!="SED" & Fraction!="PS")
psCOI_before_prune <- subset_samples(psCOI_cleaned,Fraction!="SED" & Fraction!="PS")

# investigate which samples are removed 
samples_before_prune <- colnames(otu_table(psCOI_cleaned))

# Remove samples with a read number of zero
psCOI_cleaned <- prune_samples(sample_sums(psCOI_cleaned) > 0, psCOI_cleaned)

# investigate which samples are removed 
samples_after_prune <- colnames(otu_table(psCOI_cleaned))
pruned_samples <- setdiff(samples_before_prune, samples_after_prune)

# Remove MOTUs which have a total abundance of zero after removing samples during all of the previous steps
psCOI_cleaned <- prune_taxa(rowSums(otu_table(psCOI_cleaned)) > 0, psCOI_cleaned)

# remove duplicate samples for samples retrieved in 2022 and 2023
excluded_samples <- c("DBQ_AAIDOSTA_1", "DBQ_AAIMOSTA_1", "DBQ_AAIMOSTA_2",
                      "DBQ_AAIFOSTA_2", "DBQ_AAMLOSTA_1", "DBQ_AAMNOSTA_1",
                      "DBQ_AAMPOSTA_1", "DBQ_ABBSOSTA_2", "DBQ_ABBUOSTA_2",
                      "DBQ_ABBXOSTA_2", "DBQ_ABBZOSTA_2", "DBQ_ABCBOSTA_2",
                      "DBQ_ABCDOSTA_2", "DBQ_ABGMOSTA_2", "DBQ_ABCHOSTA_2")

psCOI_cleaned  <- subset_samples(psCOI_cleaned, !(sample_names(psCOI_cleaned) %in% excluded_samples))

# Save this phyloseq object as the most unfiltered ARMS data set
saveRDS(psCOI_cleaned, file = file.path(boldigger_dir, "psCOI_unfiltered_ARMS.rds"))

# get number of remaining samples and MOTUs
psCOI_cleaned

# get number of reads
sum(sample_sums(psCOI_cleaned))

# Get percentage of MOTUs classified at phylum level
motu_phylum <- nrow(otu_table(subset_taxa(psCOI_cleaned,!is.na(phylum))))/nrow(otu_table(psCOI_cleaned))
message(round(motu_phylum*100, 2), "% of MOTUs have been classified at phylum level")

# get percentage of reads classified at phylum level
reads_phylum <- sum(sample_sums(subset_taxa(psCOI_cleaned,!is.na(phylum))))/sum(sample_sums(psCOI_cleaned))
message(round(reads_phylum*100, 2), "% of reads have been classified at phylum level")

# Get percentage of MOTUs classified at species level
motu_phylum <- nrow(otu_table(subset_taxa(psCOI_cleaned,!is.na(species))))/nrow(otu_table(psCOI_cleaned))
message(round(motu_phylum*100, 2), "% of MOTUs have been classified at species level")

# get percentage of reads classified at species level
reads_phylum <- sum(sample_sums(subset_taxa(psCOI_cleaned,!is.na(species))))/sum(sample_sums(psCOI_cleaned))
message(round(reads_phylum*100, 2), "% of reads have been classified at species level")

# create file for WoRMS taxon match service

# The unfiltered data set is the least exclusive, i.e., all MOTUs appearing in this data set are found in all potential filtered data sets
# Get a table of the unfiltered data set with MOTU names and Genus and Species assignments in one column for taxa classified down to species level

taxa <- cbind(rownames(tax_table(subset_taxa(psCOI_cleaned, !is.na(species)))), tax_table(subset_taxa(psCOI_cleaned, !is.na(species)))[, 6])
colnames(taxa) <- c("MOTU","ScientificName") # The LifeWatch e-Lab service used later on needs the Species column to be named "ScientificName"

# Write table to file which will be checked for correct species names in WoRMS via LifeWatch Belgium's e-Lab services https://www.lifewatch.be/data-services/
# https://www.marinespecies.org/msbias/aphia.php?p=match 
write.table(taxa, file.path(post_dir, "COI_species_check_WoRMS.txt"),
            sep = ";", quote = F, row.names = F)
