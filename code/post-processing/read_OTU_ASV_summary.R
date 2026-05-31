#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyverse)
library(data.table)
library(vegan)
library(ggplot2)

meta_dir <- "metadata/generated_meta"

# read phyloseq object
coi_phylo <- readRDS("output/final_taxonomy/non-rarefied_phyloseq.rds")

# read rarefied phyloseq object
rarefied_phylo <- readRDS(file = file.path("output/final_taxonomy/rarefied_phyloseq.rds"))

# read presence/absence tables and richness tables (summarized number of asvs/otus)
asv_presence_absence <- read.csv("output/stats/ASV_presence_absence.csv", 
                                sep = "\t", header = TRUE)
otu_presence_absence <- read.csv("output/stats/OTU_presence_absence.csv", 
                                sep = "\t", header = TRUE)

asv_otu <- read.csv("output/final_taxonomy/ASV_OTU_mapping.csv",
                      header = TRUE, sep = ";")

asv_count <- read.csv("output/COI/COI_ASV_counts_nosingle.txt", 
                      header = TRUE, sep = "\t")

track_miseq <- read.csv("output/COI/track_miseq.txt", header = TRUE, sep = "\t")
track_novaseq <- read.csv("output/COI/track_novaseq.txt", header = TRUE, sep = "\t")

# Read sample metadata
miseq_sample_sum <- read.table(file = file.path(meta_dir, "miseq_COI_sample_data.txt"), 
                            header = T, row.names = 1, check.names = F, sep = "\t", strip.white = T)

novaseq_sample_sum <- read.table(file = file.path(meta_dir, "novaseq_COI_sample_data.txt"), 
                            header = T, row.names = 1, check.names = F, sep = "\t", strip.white = T)

# how many OTUs were assigned at species level? 

sp <- nrow(tax_table(subset_taxa(rarefied_phylo, !is.na(species))))
otus <- nrow(tax_table(rarefied_phylo))
asvs <- dim(asv_presence_absence[, -1])

asvs

otus

sp

sp/otus

# make long format table of otus, asvs and samples (section made with help from chatgpt)
otu_long <- asv_otu %>%
  pivot_longer(
    cols = starts_with("MOTU_string"),
    names_to = "string_col",
    values_to = "ASV"
  ) %>%
  filter(!is.na(ASV))

asv_long <- asv_count %>%
  pivot_longer(
    cols = -X,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  rename(ASV = X)

merged_long <- otu_long %>%
  left_join(asv_long, by = "ASV") %>%
    filter(Abundance != 0)


# do the same sum but for year by year
years <- sort(unique(sample_data(coi_phylo)$RetrievedYear))
class(years) <- "character"


# one loop for the filtered data (before rarefaction)
for (year in years){

    # subset coi_phylo to year
    subset_phylo <- subset_samples(coi_phylo, RetrievedYear == year)
    # remove all motus with zero reads after subsetting
    otu_tab <- otu_table(subset_phylo)[apply(otu_table(subset_phylo), 1, sum) > 0, ]

    # extract MOTUs remaining and sample names
    sample_names <- sample_names(subset_phylo)
    otu_names <- rownames(otu_tab)

    # subset sample names to merged_long
    merged_subset <- merged_long %>%
        filter(Sample %in% sample_names) %>%
        filter(MOTU %in% otu_names)

    ASV_sum <- n_distinct(merged_subset$ASV)
    OTU_sum <- n_distinct(merged_subset$MOTU)

    otu_tab <- subset(otu_table(subset_phylo), rownames(otu_table(subset_phylo)) %in% rownames(otu_tab))
    otu_phylo <- merge_phyloseq(otu_tab, tax_table(subset_phylo), sample_data(subset_phylo))    

    filtered_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_tab))
    
    species_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(species))))

    message(year, ": ", ASV_sum, " ASVs, ", OTU_sum, " OTUs, ", species_sum, " species")
}

for (year in years){

    # subset coi_phylo to year
    subset_phylo <- subset_samples(rarefied_phylo, RetrievedYear == year)
    # remove all motus with zero reads after subsetting
    otu_tab <- otu_table(subset_phylo)[apply(otu_table(subset_phylo), 1, sum) > 0, ]

    # extract MOTUs remaining and sample names
    sample_names <- sample_names(subset_phylo)
    otu_names <- rownames(otu_tab)

    # subset sample names to merged_long
    merged_subset <- merged_long %>%
        filter(Sample %in% sample_names) %>%
        filter(MOTU %in% otu_names)

    ASV_sum <- n_distinct(merged_subset$ASV)
    OTU_sum <- n_distinct(merged_subset$MOTU)

    otu_tab <- subset(otu_table(subset_phylo), rownames(otu_table(subset_phylo)) %in% rownames(otu_tab))
    otu_phylo <- merge_phyloseq(otu_tab, tax_table(subset_phylo), sample_data(subset_phylo))    

    filtered_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_tab))
    
    species_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(species))))
    genera_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(genus))))
    family_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(family))))
    order_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(order))))
    class_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(class))))
    phylum_sum <- nrow(tax_table(subset_taxa(otu_phylo, !is.na(phylum))))

    message(year, ": ", ASV_sum, " ASVs, ", OTU_sum, " OTUs, ", species_sum, " species, ", genera_sum, " genera")
    message(year, ": ", family_sum, " families, ", order_sum, " orders, ", class_sum, " classes, ", phylum_sum, " phyla")
}


# check if there are any ASVs unique for sample ERR4018467 that were removed
asvs_ERR4018467 <- asv_count[asv_count$ERR4018467 == apply(asv_count[2:ncol(asv_count)], 1, sum), "X"]
number_ERR4018467 <- length(asvs_ERR4018467)

# subset the asv_otu summary to motus remaining and sum the No_ASV table
filtered_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(coi_phylo)))
sum(filtered_asv$No_ASVs)

# summarize the No_ASVs column for total number of ASVs before and after rarefaction
rarefied_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(rarefied_phylo)))
# remove the number of asvs unique for sample ERR4018467 for the rarefied sum 
sum(rarefied_asv$No_ASVs) - number_ERR4018467




######## calculate nr of input and output sequences #########
miseq_input <- sum(track_miseq$input)
miseq_input

novaseq_input <- sum(track_novaseq$input)
novaseq_input

total_input <- miseq_input + novaseq_input
total_input

########## calculate number of reads, OTUs, species remaining before rarefaction ##########
coi_phylo
reads_sum <- sum(sample_sums(coi_phylo))
reads_sum 
# otus and species before rarefaction
dim(otu_table(coi_phylo))
dim(otu_table(subset_taxa(coi_phylo, !is.na(species))))

miseq_bf <- subset_samples(coi_phylo, SeqMethod == "miseq")
miseq_bf <- prune_taxa(taxa_sums(miseq_bf) > 0, miseq_bf)
novaseq_sum <- sum(sample_sums(miseq_bf))
dim(otu_table(miseq_bf))
dim(otu_table(subset_taxa(miseq_bf, !is.na(species))))

novaseq_bf <- subset_samples(coi_phylo, SeqMethod == "novaseq")
novaseq_bf <- prune_taxa(taxa_sums(novaseq_bf) > 0, novaseq_bf)
novaseq_sum <- sum(sample_sums(novaseq_bf))
novaseq_sum
dim(otu_table(novaseq_bf))
dim(otu_table(subset_taxa(novaseq_bf, !is.na(species))))

########## calculate number of reads, OTUs, species remaining AFTER rarefaction ##########

# calculate number of reads remaining after rarefaction
reads_sum <- sum(sample_sums(rarefied_phylo))
reads_sum
# otus & species after rarefaction
dim(otu_table(rarefied_phylo))
dim(otu_table(subset_taxa(rarefied_phylo, !is.na(species))))

miseq_af <- subset_samples(rarefied_phylo, SeqMethod == "miseq")
miseq_af <- prune_taxa(taxa_sums(miseq_af) > 0, miseq_af)
miseq_sum <- sum(sample_sums(miseq_af))
miseq_sum
dim(otu_table(miseq_af))
dim(otu_table(subset_taxa(miseq_af, !is.na(species))))

novaseq_af <- subset_samples(rarefied_phylo, SeqMethod == "novaseq")
novaseq_af <- prune_taxa(taxa_sums(novaseq_af) > 0, novaseq_af)
novaseq_sum <- sum(sample_sums(novaseq_af))
novaseq_sum
dim(otu_table(novaseq_af))
dim(otu_table(subset_taxa(novaseq_af, !is.na(species))))


############# calculate no asvs before and after rarefaction ##############

# asv miseq before rarefaction 
filtered_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(miseq_bf)))
sum(filtered_asv$No_ASVs)

# asv miseq after rarefaction 
rarefied_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(miseq_af)))
# remove the number of asvs unique for sample ERR4018467 for the rarefied sum 
sum(rarefied_asv$No_ASVs) - number_ERR4018467

# asv novaseq before rarefaction 
filtered_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(novaseq_bf)))
sum(filtered_asv$No_ASVs)

# asv novaseq after rarefaction 
rarefied_asv <- subset(asv_otu, asv_otu$MOTU %in% rownames(otu_table(novaseq_af)))
# remove the number of asvs unique for sample ERR4018467 for the rarefied sum 
sum(rarefied_asv$No_ASVs)

################# calculate number of reads year by year bf and af raref. #################
bf_2019 <- subset_samples(coi_phylo, RetrievedYear == "2019")
bf_2019 <- prune_taxa(taxa_sums(bf_2019) > 0, bf_2019)
sum(sample_sums(bf_2019))

af_2019 <- subset_samples(rarefied_phylo, RetrievedYear == "2019")
af_2019 <- prune_taxa(taxa_sums(af_2019) > 0, af_2019)
sum(sample_sums(af_2019))

bf_2020 <- subset_samples(coi_phylo, RetrievedYear == "2020")
bf_2020 <- prune_taxa(taxa_sums(bf_2020) > 0, bf_2020)
sum(sample_sums(bf_2020))

af_2020 <- subset_samples(rarefied_phylo, RetrievedYear == "2020")
af_2020 <- prune_taxa(taxa_sums(af_2020) > 0, af_2020)
sum(sample_sums(af_2020))

bf_2021 <- subset_samples(coi_phylo, RetrievedYear == "2021")
bf_2021 <- prune_taxa(taxa_sums(bf_2021) > 0, bf_2021)
sum(sample_sums(bf_2021))

af_2021 <- subset_samples(rarefied_phylo, RetrievedYear == "2021")
af_2021 <- prune_taxa(taxa_sums(af_2021) > 0, af_2021)
sum(sample_sums(af_2021))

bf_2022 <- subset_samples(coi_phylo, RetrievedYear == "2022")
bf_2022 <- prune_taxa(taxa_sums(bf_2022) > 0, bf_2022)
sum(sample_sums(bf_2022))

af_2022 <- subset_samples(rarefied_phylo, RetrievedYear == "2022")
af_2022 <- prune_taxa(taxa_sums(af_2022) > 0, af_2022)
sum(sample_sums(af_2022))

bf_2023 <- subset_samples(coi_phylo, RetrievedYear == "2023")
bf_2023 <- prune_taxa(taxa_sums(bf_2023) > 0, bf_2023)
sum(sample_sums(bf_2023))

af_2023 <- subset_samples(rarefied_phylo, RetrievedYear == "2023")
af_2023 <- prune_taxa(taxa_sums(af_2023) > 0, af_2023)
sum(sample_sums(af_2023))

# calculate number of reads from species level assignment
reads_sum_sp <- sum(sample_sums(subset_taxa(rarefied_phylo, !is.na(species))))
reads_sum_sp

reads_sum_sp/reads_sum


# how many % of reads did the 183 removed OTUs represent?
removed_otus <- setdiff(rownames(otu_table(coi_phylo)), rownames(otu_table(rarefied_phylo)))

sum(sample_sums(prune_taxa(removed_otus, coi_phylo)))/sum(sample_sums(coi_phylo))*100

# did any species get removed from rarefaction? 
nrow(tax_table(subset_taxa(coi_phylo, !is.na(species)))) - nrow(tax_table(subset_taxa(rarefied_phylo, !is.na(species))))

1 - nrow(tax_table(subset_taxa(rarefied_phylo, !is.na(species))))/nrow(tax_table(subset_taxa(coi_phylo, !is.na(species))))


############# summarize read numbers, nOTU and nASV per sample ##############
trans_otu <- as.data.frame(t(otu_presence_absence))
colnames(trans_otu) <- trans_otu[1, ]
trans_otu <- trans_otu[-1, ]

perma_df_otu <- merge(as.data.frame(sample_data(rarefied_phylo)), trans_otu, by = "row.names")
rownames(perma_df_otu) <- perma_df_otu[, 1]
perma_df_otu <- perma_df_otu[, -1]

trans_asv <- as.data.frame(t(asv_presence_absence))
colnames(trans_asv) <- trans_asv[1, ]
trans_asv <- trans_asv[-1, ]

perma_df_asv <- merge(as.data.frame(sample_data(rarefied_phylo)), trans_asv, by = "row.names")
rownames(perma_df_asv) <- perma_df_asv[, 1]
perma_df_asv <- perma_df_asv[, -1]

data.frame(sample_sums(coi_phylo), sample_data(coi_phylo)[, 1])
data.frame(sample_sums(rarefied_phylo), sample_data(rarefied_phylo)[, 1])

perma_df_otu <- perma_df_otu[order(perma_df_otu$SeqMethod), ]
otu_matrix <- perma_df_otu[, 10:ncol(perma_df_otu)]
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric"
data.frame(rowSums(otu_matrix), perma_df_otu$MaterialSampleID)

perma_df_asv <- perma_df_asv[order(perma_df_asv$SeqMethod), ]
asv_matrix <- perma_df_asv[, 10:ncol(perma_df_asv)]
asv_matrix <- as.matrix(asv_matrix)
storage.mode(asv_matrix) <- "numeric"
data.frame(rowSums(asv_matrix), perma_df_asv$MaterialSampleID)


#################### track reads year by year ##################

samples19 <- rownames(miseq_sample_sum[miseq_sample_sum$DeployedDate == "180418", ])
reads19 <- colSums(track_miseq[track_miseq$X %in% samples19, -1])
reads19/reads19[1]

samples20 <- rownames(miseq_sample_sum[miseq_sample_sum$RetrievedYear == "2020", ])
reads20 <- colSums(track_miseq[track_miseq$X %in% samples20, -1])
reads20/reads20[1]

samples21 <- rownames(miseq_sample_sum[miseq_sample_sum$RetrievedYear == "2021", ])
reads21 <- colSums(track_miseq[track_miseq$X %in% samples21, -1])
reads21/reads21[1]

samples22 <- rownames(novaseq_sample_sum[novaseq_sample_sum$RetrievedYear == "2022", ])
reads22 <- colSums(track_novaseq[track_novaseq$X %in% samples22, -1])
reads22/reads22[1]

samples23 <- rownames(novaseq_sample_sum[novaseq_sample_sum$RetrievedYear == "2023", ])
reads23 <- colSums(track_novaseq[track_novaseq$X %in% samples23, -1])
reads23/reads23[1]

reads19+reads20+reads21+reads22+reads23

