#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(ggplot2)
library(tibble)
library(scales)

# specify path to figures output
figures_dir <- "results/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)
nmds_dir <- "results/figures/nmds_plots"
if(!dir.exists(nmds_dir)) dir.create(nmds_dir)
scatter_dir <- "results/figures/scatter_plots"
if(!dir.exists(scatter_dir)) dir.create(scatter_dir)

nmds_ambi_dir <- "results/figures/nmds_plots/ambi"
if(!dir.exists(nmds_ambi_dir)) dir.create(nmds_ambi_dir)
scatter_ambi_dir <- "results/figures/scatter_plots/ambi"
if(!dir.exists(scatter_ambi_dir)) dir.create(scatter_ambi_dir)

# specify path to richness plots
richness_dir <- file.path(figures_dir, "richness_plots")
if(!dir.exists(richness_dir)) dir.create(richness_dir)

# specify path to final taxonomy dir
post_dir <- "output/final_taxonomy"

# read final taxonomy rds
# read phyloseq object
rarefied_phylo <- readRDS(file = file.path("output/final_taxonomy/rarefied_phyloseq.rds"))

# load species list from sample site, alien species list and red list
skagerrak_table <- read.csv("output/final_taxonomy/TaxaExport_2026-04-02_16.32.22.csv",
                      header = TRUE, sep = ";")
alien_list <- read.csv("output/final_taxonomy/TaxaExport_2026-04-07_11.12.36_Frammande_i_Sverige.csv",
                        header = TRUE, sep = ";")
red_list <- read.csv("output/final_taxonomy/TaxaExport_2026-04-07_11.15.51_RE_CR_EN_VU.csv",
                        header = TRUE, sep = ";")
# load ambi very sensitive to disturbance list from lifewatch
traits_table <- read.csv("output/final_taxonomy/LW_Traits_data___downloaded_on_2026-04-07-15-58.tab",
                      header = TRUE, sep = "\t")
# load ASV OTU mapping csv
asv_otu <- read.csv("output/final_taxonomy/ASV_OTU_mapping.csv",
                      header = TRUE, sep = ";")


######### calculate mean ASV and OTU count for each chosen category and plot ##############

################# create AMBI very sensitive to disturbance subset ###############
traits_species <- paste(traits_table$Genus, traits_table$Species)
traits_common <- intersect(tax_table(rarefied_phylo)[, 6], traits_species)
traits_phylo <- subset_taxa(rarefied_phylo, species %in% traits_common)

# subset and sort asv_otu table to MOTU and No_ASVs
asv_subset <- data.frame(asv_otu[, 3])
colnames(asv_subset) <- "No_ASVs"
rownames(asv_subset) <- asv_otu[, 1]

# merge otu_table and asv_subset
asvtab_subset <- subset(asv_subset, rownames(asv_subset) %in% rownames(otu_table(traits_phylo)))
merge_otu_asv <- merge(asvtab_subset, otu_table(traits_phylo), by = 0)

logical_subset <- merge_otu_asv[, -c(1, 2)] != 0
rownames(logical_subset) <- merge_otu_asv[, 1]

# create presence/absence table of OTUs
presence_absence_otus <- 1 * logical_subset

# produce a data set with the number of ASVs in each (non-zero) OTU for each sample
for (i in 1:nrow(logical_subset)){
  logical_subset[i, ] <- ifelse(logical_subset[i, ] == 1, merge_otu_asv$No_ASVs[i], 0)
}


################# make mean OTU/ASV plots #########################

# obtain all sub-categories of chosen category
meta_cats <- sort(unique(sample_data(traits_phylo)$RetrievedYear))
class(meta_cats) <- "character"

# create empty data base
sum_summary <- data.frame(
  OTU_sum = c(),
  ASV_sum = c()
)
mean_summary <- data.frame(
  OTU_mean = c(),
  OTU_sd = c(),
  ASV_mean = c(),
  ASV_sd = c()
)

for (cat in meta_cats){
  current_cat <- cat
  # create phyloseq object of the specific sub-category. e.g. 2019 or 2020
  cat_subset <- subset_samples(traits_phylo, RetrievedYear == current_cat)
  # extract sample names from samples present in the subset
  cat_samples <- sample_names(cat_subset)

  # subset the samples to the ASV/OTU logical subset
  cat_logic <- logical_subset[, cat_samples]

  # calculate mean OTUs per sample within the sub-category 
  OTU_sum <- colSums(cat_logic != 0, na.rm = TRUE)
  OTU_mean <- mean(OTU_sum, na.rm = TRUE)
  OTU_sd <- sd(OTU_sum)

  # calculate the mean ASV number per sample within the sub-category
  ASV_sum <- colSums(cat_logic, na.rm = TRUE)
  ASV_mean <- mean(ASV_sum, na.rm = TRUE)
  ASV_sd <- sd(ASV_sum)

  # collect ASV_sum, OTU_sum, ASV_mean and OTU_mean in a data base
    output1 <- data.frame(
      OTU_sum = OTU_sum,
      ASV_sum = ASV_sum)
    output2 <- data.frame(
      OTU_mean = OTU_mean,
      OTU_sd = OTU_sd,
      ASV_mean = ASV_mean, 
      ASV_sd = ASV_sd, 
      row.names = as.character(cat))
  sum_summary <- rbind(sum_summary, output1)
  mean_summary <- rbind(mean_summary, output2)
}

mean_summary <- rownames_to_column(mean_summary, "Year") 

# investigate different plots (barplots, regular plots with trendlines, ggplot2)

# mean OTU plot with standard deviation as bars

# plot the point plot
OTU_mean_plot <- ggplot(mean_summary, aes(x = Year, y = OTU_mean)) + 
  geom_point()+
  geom_errorbar(aes(ymin = OTU_mean - OTU_sd, ymax = OTU_mean + OTU_sd), width = .1,
                position = position_dodge(0.05), color = "#1e1e21") +
  theme_minimal() + 
      theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank()
  )

ggsave(plot = OTU_mean_plot, file = file.path(scatter_ambi_dir, "AMBI_mean_sd_OTUs.png"),
       dpi = 400)

  # plot the point plot
ASV_mean_plot <- ggplot(mean_summary, aes(x = Year, y = ASV_mean)) + 
  geom_point()+
  geom_errorbar(aes(ymin = ASV_mean - ASV_sd, ymax = ASV_mean + ASV_sd), width = .1,
                position = position_dodge(0.05), color = "#1e1e21") +
  theme_minimal() + 
      theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank()
  )

  ggsave(plot = ASV_mean_plot, file = file.path(scatter_ambi_dir, "AMBI_mean_sd_ASVs.png"),
       dpi = 400)



################# make NMDS (OTU) plots #########################

nmds_presence_absence <- t(presence_absence_otus)
nmds_presence_absence <- as.matrix(nmds_presence_absence)

nmds <- metaMDS(nmds_presence_absence, distance = "jaccard")
# stressplot(nmds)

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
data_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
data_scores <- transform(merge(sample_data(traits_phylo)[, 3], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(traits_phylo)[, 5], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(traits_phylo)[, 8], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(traits_phylo)[, 9], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
colors <- c("#54accb", "#723ed2", "#f353ac", "#fdcb31")
gradient_colors <- colorRampPalette(colors)(5)
# show_col(gradient_colors)


# plot nmds plot with RetrievedYear and SeqMethod
nmds_year_seq <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "SeqMethod")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_year_seq, file = file.path(nmds_ambi_dir, "AMBI_NMDS_OTU_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_year_unit <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "ARMSUnit")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_year_unit, file = file.path(nmds_ambi_dir, "AMBI_NMDS_OTU_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_year_frac <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "FractionGroup")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_year_frac, file = file.path(nmds_ambi_dir, "AMBI_NMDS_OTU_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



################# make NMDS (ASV) plots #########################

nmds_asv <- t(logical_subset)

nmds <- metaMDS(nmds_asv, distance = "bray")
# stressplot(nmds)
# nmds

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
asv_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
asv_scores <- transform(merge(sample_data(traits_phylo)[, 3], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(traits_phylo)[, 5], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(traits_phylo)[, 8], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(traits_phylo)[, 9], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
colors <- c("#54accb", "#723ed2", "#f353ac", "#fdcb31")
gradient_colors <- colorRampPalette(colors)(5)
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
nmds_asv_year_seq <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "SeqMethod")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_asv_year_seq, file = file.path(nmds_ambi_dir, "AMBI_NMDS_ASV_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_asv_year_unit <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "ARMSUnit")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_asv_year_unit, file = file.path(nmds_ambi_dir, "AMBI_NMDS_ASV_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_asv_year_frac <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  scale_fill_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "RetrievedYear", y = "NMDS2", shape = "FractionGroup")  + 
  theme_minimal() + 
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
    theme(
      plot.background = element_rect(fill = "#fdfbf7"),
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
      panel.grid.minor = element_blank())

ggsave(plot = nmds_asv_year_frac, file = file.path(nmds_ambi_dir, "AMBI_NMDS_ASV_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



