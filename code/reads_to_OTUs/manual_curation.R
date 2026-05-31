#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyverse)
library(data.table)
library(vegan)
library(ggplot2)

# reference and create figures dir
# specify path to figures output
figures_dir <- "results/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)
# specify path to richness plots
richness_dir <- file.path(figures_dir, "richness_plots")
if(!dir.exists(richness_dir)) dir.create(richness_dir)

# specify path to output dir
post_dir <- "output/final_taxonomy"

# read phyloseq object
coi_phylo <- readRDS("output/COI/BOLDigger/psCOI_unfiltered_ARMS.rds")

# load worms taxon matches 
worms_match <- read.csv("output/final_taxonomy/coi_species_check_worms_matched.txt",
                      header = TRUE, sep = ";")

# load species list from Artdatabanken
skagerrak_table <- read.csv("output/final_taxonomy/TaxaExport_2026-04-02_16.32.22.csv",
                      header = TRUE, sep = ";")


########## change the MaterialSampleID column to format Year_ARMSUnit_FractionGroup ##########
year <- sample_data(coi_phylo)[, 8]
unit <- sample_data(coi_phylo)[, 9]
fraction <- sample_data(coi_phylo)[, 5]

material_sample <- paste(year$RetrievedYear, unit$ARMSUnit, fraction$FractionGroup, sep = "_")
sample_data(coi_phylo)[, 1] <- material_sample

########## remove contaminants, insects and spiders ##########
# investigated species list and identified contaiminant taxa
species_phylo <- subset_taxa(coi_phylo, !is.na(species)) 
# View(tax_table(species_phylo))
# unique(tax_table(coi_phylo)[, 1])

coi_phylo  <- subset_taxa(coi_phylo, class != "Insecta")
coi_phylo  <- subset_taxa(coi_phylo, class != "Arachnida")
contaminant_phyla <- c("Pseudomonadota", "Basidiomycota", "Ascomycota", "Mucoromycota")
coi_phylo  <- subset_taxa(coi_phylo, !(phylum %in% contaminant_phyla))
contaminant_species <- c("Homo sapiens", "Sus scrofa")
coi_phylo  <- subset_taxa(coi_phylo, !(species %in% contaminant_species))

######### save final filtered phyloseq object before rarefaction #########
saveRDS(coi_phylo, file = file.path(post_dir, "non-rarefied_phyloseq.rds"))

########## rarefy to 6 761 reads #########

# visualize rarefaction to 10 000 reads 
tab <- otu_table(coi_phylo)
class(tab) <- "matrix" # as.matrix() will do nothing
tab <- t(tab)
png(file.path(figures_dir, "rarefaction/10000_reads_rarefaction.png"), res = 400, width = 16, height = 12, units = "in")
rarecurve(tab, step = 500, sample = 10000)
title(main = "Example rarefaction curves")
dev.off()

# visualiza rarefaction to 8 000 reads 
tab <- otu_table(coi_phylo)
class(tab) <- "matrix" # as.matrix() will do nothing
tab <- t(tab)
png(file.path(figures_dir, "rarefaction/8000_reads_rarefaction.png"), res = 400, width = 16, height = 12, units = "in")
rarecurve(tab, step = 500, sample = 8000)
title(main = "Example rarefaction curves")
dev.off()

tab <- otu_table(coi_phylo)
class(tab) <- "matrix" # as.matrix() will do nothing
tab <- t(tab)
png(file.path(figures_dir, "rarefaction/6761_reads_rarefaction.png"), res = 400, width = 16, height = 12, units = "in")
rarecurve(tab, step = 500, sample = 6761, label = TRUE, ylab = "OTUs", xlab = "Sequencing depth")
dev.off()

coi_rarecurve <- rarecurve(tab, step = 500, sample = 6761, tidy = TRUE, label = TRUE, ylab = "OTUs", xlab = "Sequencing depth")

prune_phylo <- prune_samples(!(sample_names(species_phylo) %in% c("DBQ_ABBTOSTA_2", "DBQ_AAMOSTA_1", "DBQ_AAMKOSTA_1",
                                                 "DBQ_ABBWOSTA_2", "DBQ_ABCGOSTA_2", "DBQ_AAMOOSTA_1", "DBQ_ABCAOSTA_2",
                                                 "DBQ_ABCCOSTA_2", "DBQ_AAIGOSTA_2", "DBQ_AAMMOSTA_1", "DBQ_ABBQOSTA_2",
                                                 "DBQ_AAICOSTA_1", "DBQ_ABCEOSTA_2")), species_phylo)
# visualize rarefaction to 6 761 reads 
tab <- otu_table(prune_phylo)
class(tab) <- "matrix" # as.matrix() will do nothing
tab <- t(tab)
png(file.path(figures_dir, "rarefaction/6761_reads_rarefaction_pruned.png"), res = 400, width = 16, height = 12, units = "in")
rarecurve(tab, step = 500, sample = 6761, label = TRUE, ylab = "OTUs", xlab = "Sequencing depth")
dev.off()

# rarefy to 6 761 reads
rarefied_phylo <- rarefy_even_depth(coi_phylo, rngseed = 1, sample.size = 6761, replace = F)

########### taxon match service ###########

# match remaining phyloseq species rows to input ScientificName column in taxon match service file
# this is to only include the species names present in our filtered data set

# In some cases, there were some fuzzy or non-exact matches of our taxa names vs. the ones found in WoRMS
# In such cases, the entry in the column "accepted_name_aphia_worms" will be empty
# replace these entries with our original names from the "scientificname" column
worms_subset <- worms_match[, c("ScientificName", "ScientificName_accepted")]

worms_subset$ScientificName_accepted <- ifelse(worms_subset$ScientificName_accepted == "", worms_subset$ScientificName, worms_subset$ScientificName_accepted)

# if there is a bracket in the ScientificName_accepted column, replace that name with ScientificName
worms_subset$ScientificName_accepted <- ifelse(grepl("\\(", worms_subset$ScientificName_accepted), worms_subset$ScientificName, worms_subset$ScientificName_accepted)

# if the ScientificName_accepted has three characters (denotes sub-species), change to ScientificName
worms_subset$ScientificName_accepted <- ifelse(lengths(strsplit(worms_subset$ScientificName_accepted, " ")) > 2, worms_subset$ScientificName, worms_subset$ScientificName_accepted)

subset_taxa(rarefied_phylo, !is.na(species)) 
 
for (species in tax_table(rarefied_phylo)[, 6]){
  species_row <- grep(species, tax_table(rarefied_phylo)[, 6])
  tax_table(rarefied_phylo)[species_row, 6] <- ifelse(length(grep(species, worms_subset$ScientificName_accepted)) == 0, species, grep(species, worms_subset$ScientificName_accepted, value = TRUE))
}

subset_taxa(rarefied_phylo, !is.na(species)) 

species_reads <- cbind(tax_table(rarefied_phylo)[, 6], otu_table(rarefied_phylo))

write.table(species_reads, file = file.path(post_dir, "filtered_species_list.csv"),
            quote = FALSE, sep = ";", row.names = FALSE)

# save rarefied and taxon match data set to rds
saveRDS(rarefied_phylo, file = file.path(post_dir, "rarefied_phyloseq.rds"))


######### plot rarefaction curves using ggplot2 ##########

meta <- data.frame(sample_data(coi_phylo), stringsAsFactors = FALSE)
meta$Site <- sample_names(coi_phylo)

plot_df <- coi_rarecurve %>% 
  left_join(meta, by = "Site")

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
rare_plot <- ggplot(plot_df, aes(x = Sample, y = Species, color = RetrievedYear, group = MaterialSampleID)) + 
  scale_color_manual(values = gradient_colors) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = 6761, color = "#1e1e21", size = 0.3) +
  xlim(0, 50000) +
  facet_grid(~RetrievedYear) +
  labs(x = "Sequencing depth", colour = "Retrieved Year", y = "Number of OTUs")  + 
  theme_minimal() + 
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),

      strip.text = element_text(size = 20, color = "#1e1e21"),

      legend.text = element_text(size = 14, color = "#1e1e21"),  
      legend.title = element_text(size = 18, color = "#1e1e21"),

      legend.box.just = "center",
      legend.justification = "center",
        
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 20),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
      axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank())


ggsave(plot = rare_plot, file = file.path(figures_dir, "rarefaction/6761_reads_rarefaction.png"),
       dpi = 400, height = 10.5, width = 14)

