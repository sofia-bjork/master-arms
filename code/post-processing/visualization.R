#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(ggplot2)
library(tibble)
library(scales)

miseq <- sum_summary[sum_summary$seq_method == "miseq", ]
novaseq <- sum_summary[sum_summary$seq_method == "novaseq", ]

set.seed(1)  

# specify path to figures output
figures_dir <- "results/figures"
if(!dir.exists(figures_dir)) dir.create(figures_dir)
nmds_dir <- "results/figures/nmds_plots"
if(!dir.exists(nmds_dir)) dir.create(nmds_dir)
scatter_dir <- "results/figures/scatter_plots"
if(!dir.exists(scatter_dir)) dir.create(scatter_dir)

# specify path to richness plots
richness_dir <- file.path(figures_dir, "richness_plots")
if(!dir.exists(richness_dir)) dir.create(richness_dir)

# specify path to final taxonomy dir
post_dir <- "output/final_taxonomy"

# specify path to presence/absence and richness directory 
stats_dir <- "output/stats"
if(!dir.exists(stats_dir)) dir.create(stats_dir)

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
# load ASV counts table from after chimera removal
asv_count <- read.csv("output/COI/COI_ASV_counts_nosingle.txt", 
                      header = TRUE, sep = "\t")

######### calculate mean ASV and OTU count for each chosen category and plot ##############

# subset and sort asv_otu table to MOTU and No_ASVs
asv_subset <- data.frame(asv_otu[, 3])
colnames(asv_subset) <- "No_ASVs"
rownames(asv_subset) <- asv_otu[, 1]

# merge otu_table and asv_subset
asvtab_subset <- subset(asv_subset, rownames(asv_subset) %in% rownames(otu_table(rarefied_phylo)))
merge_otu_asv <- merge(asvtab_subset, otu_table(rarefied_phylo), by = 0)

logical_subset <- merge_otu_asv[, -c(1, 2)] != 0
rownames(logical_subset) <- merge_otu_asv[, 1]

# create presence/absence table of OTUs
presence_absence_otus <- 1 * logical_subset

# produce a logical data set with number of asvs present in the samples
for (i in 1:nrow(logical_subset)){
  logical_subset[i, ] <- ifelse(logical_subset[i, ] == 1, merge_otu_asv$No_ASVs[i], 0)
}

# collect ASVs still present in the data after otu clustering, lulu curation , taxonomic assignment and firther filtering
all_asvs <- c()
for (row in 1:nrow(asv_otu)){
  row_asvs <- grep("ASV", asv_otu[row, 4:ncol(asv_otu)], value = TRUE)
  row_asvs <- as.character(row_asvs)
  all_asvs <- c(all_asvs, row_asvs)
}

all_asvs <- unique(all_asvs)

# filter the imported asv table based on remaining asvs
asv_presence_absence <- asv_count[asv_count$X %in% all_asvs, ]
rownames(asv_presence_absence) <- asv_presence_absence[, 1]
asv_presence_absence <- asv_presence_absence[, -1]
asv_presence_absence[asv_presence_absence > 0] <- 1

# filter the imported asv table based on samples remaning
asv_presence_absence <- asv_presence_absence[, colnames(asv_presence_absence) %in% sample_names(rarefied_phylo)]

# remove asvs and otus with zero reads remaining
asv_presence_absence <- asv_presence_absence[rowSums(asv_presence_absence[, -1]) > 0, ]
presence_absence_otus <- presence_absence_otus[rowSums(presence_absence_otus[, -1]) > 0, ]

# we now have two presence/absence tables, one for OTUs and one for ASVs
# asv_presence_absence and presence_absence_otus
# summarize them for the number of ASVs/OTUs in each sample
asv_richness <- apply(asv_presence_absence, 2, sum)
otu_richness <- apply(presence_absence_otus, 2, sum)

# save otu and asv presence/absence tables and richness tables
write.table(asv_presence_absence, file = file.path(stats_dir, "ASV_presence_absence.csv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
write.table(presence_absence_otus, file = file.path(stats_dir, "OTU_presence_absence.csv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
write.table(asv_richness, file = file.path(stats_dir, "ASV_richness.csv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
write.table(otu_richness, file = file.path(stats_dir, "OTU_richness.csv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

################# make mean OTU/ASV plots #########################

# obtain all sub-categories of chosen category
meta_cats <- sort(unique(sample_data(rarefied_phylo)$RetrievedYear))
class(meta_cats) <- "character"

otu_richness_df <- data.frame(names(otu_richness), otu_richness)

miseq <- grep("ERR", otu_richness_df$names.otu_richness., value = TRUE)
novaseq <- grep("DBQ", otu_richness_df$names.otu_richness., value = TRUE)

otu_richness_df <- otu_richness_df[rownames(otu_richness_df) %in% novaseq, ]
otu_richness_df <- otu_richness_df[rownames(otu_richness_df) %in% miseq, ]
otu_richness_df <- data.frame(names(otu_richness), otu_richness)

# Basic Density Plot
#otu_richness_df %>%
#  ggplot(aes(x = otu_richness)) +
#  # Map histogram y-axis to density scale
#  geom_histogram(aes(y = after_stat(density)), 
#                 bins = 30,
#                 fill = "lightblue",
#                 color = "white") +
#  geom_density(color = "red", linewidth = 1) +
#  theme_test() +
#  labs(title = "Density + Histogram Combination",
#       subtitle = "Overlaying smooth trends on raw data counts",
#       x = "Sepal Length (cm)")

# create empty data base
sum_summary <- data.frame(
  Year = c(),
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
  cat_subset <- subset_samples(rarefied_phylo, RetrievedYear == current_cat)
  # extract sample names from samples present in the subset
  cat_samples <- sample_names(cat_subset)

  # subset the samples to the ASV/OTU sums
  OTU_sum <- otu_richness[cat_samples]
  ASV_sum <- asv_richness[cat_samples]

  # calculate mean OTUs per sample within the sub-category 
  OTU_mean <- mean(OTU_sum, na.rm = TRUE)
  OTU_sd <- sd(OTU_sum)

  # calculate the mean ASV number per sample within the sub-category
  ASV_mean <- mean(ASV_sum, na.rm = TRUE)
  ASV_sd <- sd(ASV_sum)

  year <- rep(current_cat, length(OTU_sum))
  names(year) <- names(OTU_sum)

  # collect ASV_sum, OTU_sum, ASV_mean and OTU_mean in a data base
    output1 <- data.frame(
      Year = year,
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

seq_method <- ifelse(rownames(sum_summary) %in% miseq, "miseq", "novaseq")
sum_summary <- cbind(sum_summary, seq_method)

# investigate different plots (barplots, regular plots with trendlines, ggplot2)

# mean OTU plot with standard deviation as bars

# plot the point plot
OTU_mean_plot <- ggplot(sum_summary, aes(x = Year, y = OTU_sum, fill = seq_method)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
  labs(x = "Year", y = "OTU richness", fill = "Sequencing method") +
  scale_fill_manual(values = c("#9a77b5", "#fc8e6d")) +
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


ggsave(plot = OTU_mean_plot, file = file.path(scatter_dir, "mean_sd_OTUs.png"),
       dpi = 400)

# plot boxplot between miseq and novaseq (OTU)

# plot the point plot
seqmethod_mean_plot <- ggplot(sum_summary, aes(x = seq_method, y = OTU_sum, fill = seq_method)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
  labs(x = "Sequencing method", y = "OTU richness") +
  scale_fill_manual(values = c("#9a77b5", "#fc8e6d")) +
  theme_minimal() + 
      theme(
      legend.position = "none", 
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

ggsave(plot = seqmethod_mean_plot, file = file.path(scatter_dir, "SeqMethod_mean_sd_OTUs.png"),
       dpi = 400, width = 7, height = 7)


  # plot the boxplot
ASV_mean_plot <- ggplot(sum_summary, aes(x = Year, y = ASV_sum, fill = seq_method)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
  labs(x = "Year", y = "ASV richness", fill = "Sequencing method") +
  scale_fill_manual(values = c("#9a77b5", "#fc8e6d")) +
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

  ggsave(plot = ASV_mean_plot, file = file.path(scatter_dir, "mean_sd_ASVs.png"),
       dpi = 400)


# plot boxplot between miseq and novaseq (ASV)

# plot the point plot
asv_seqmethod_mean_plot <- ggplot(sum_summary, aes(x = seq_method, y = ASV_sum, fill = seq_method)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
  labs(x = "Sequencing method", y = "ASV richness") +
  scale_fill_manual(values = c("#9a77b5", "#fc8e6d")) +
  theme_minimal() + 
      theme(
      legend.position = "none", 
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

ggsave(plot = asv_seqmethod_mean_plot, file = file.path(scatter_dir, "SeqMethod_mean_sd_ASVs.png"),
       dpi = 400, width = 7, height = 7)



################# make NMDS (OTU) plots #########################

nmds_presence_absence <- t(presence_absence_otus)
nmds_presence_absence <- as.matrix(nmds_presence_absence)

nmds <- metaMDS(nmds_presence_absence, distance = "jaccard")
# stressplot(nmds)
nmds_stress <- nmds$stress

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
data_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
data_scores <- transform(merge(sample_data(rarefied_phylo)[, 3], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(rarefied_phylo)[, 5], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(rarefied_phylo)[, 8], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(rarefied_phylo)[, 9], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(rarefied_phylo)[, 2], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)
gradient_colors <- c("#33a8c7",
  "#52e3e1",
  "#a0e426",
  "#fdf148",
  "#ffab00",
  "#f77976",
  "#f050ae",
  "#d883ff",
  "#024fa1")

# plot nmds plot with Sequenced and SeqMethod
nmds_year_sequenced <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = RetrievedYear, colour = Sequenced)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "Sequenced", y = "NMDS2", shape = "RetrievedYear")  + 
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

ggsave(plot = nmds_year_sequenced, file = file.path(nmds_dir, "NMDS_OTU_RetrievedYear_Sequenced.png"),
       dpi = 400, width  = 8.76, height = 6.8)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)


# plot nmds plot with RetrievedYear and SeqMethod
nmds_year_seq <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_year_seq, file = file.path(nmds_dir, "NMDS_OTU_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_year_unit <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_year_unit, file = file.path(nmds_dir, "NMDS_OTU_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_year_frac <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_year_frac, file = file.path(nmds_dir, "NMDS_OTU_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



################# make NMDS (ASV) plots #########################

nmds_asv <- t(asv_presence_absence)

nmds <- metaMDS(nmds_asv, distance = "jaccard")
# stressplot(nmds)
nmds_stress <- nmds$stress

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
asv_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
asv_scores <- transform(merge(sample_data(rarefied_phylo)[, 3], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(rarefied_phylo)[, 5], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(rarefied_phylo)[, 8], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(rarefied_phylo)[, 9], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(rarefied_phylo)[, 2], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

gradient_colors <- c("#33a8c7",
  "#52e3e1",
  "#a0e426",
  "#fdf148",
  "#ffab00",
  "#f77976",
  "#f050ae",
  "#d883ff",
  "#024fa1")

# plot nmds plot with Sequenced and SeqMethod
nmds_asv_year_sequenced <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = RetrievedYear, colour = Sequenced)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
  labs(x = "NMDS1", colour = "Sequenced", y = "NMDS2", shape = "RetrievedYear")  + 
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

ggsave(plot = nmds_asv_year_sequenced, file = file.path(nmds_dir, "NMDS_ASV_RetrievedYear_Sequenced.png"),
       dpi = 400, width  = 8.76, height = 6.8)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
nmds_asv_year_seq <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_asv_year_seq, file = file.path(nmds_dir, "NMDS_ASV_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_asv_year_unit <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_asv_year_unit, file = file.path(nmds_dir, "NMDS_ASV_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
nmds_asv_year_frac <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = nmds_asv_year_frac, file = file.path(nmds_dir, "NMDS_ASV_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



################# do the same thing on a species subset ################
species_subset <- subset_taxa(rarefied_phylo, !(is.na(species)))

# merge otu_table and asv_subset
sp_otu_subset <- subset(presence_absence_otus, rownames(presence_absence_otus) %in% rownames(otu_table(species_subset)))

# make MOTU to species translation
species_MOTU <- as.data.frame(tax_table(species_subset)[, 6])
logical_species <- merge(species_MOTU, sp_otu_subset, by = "row.names")
rownames(logical_species) <- logical_species[, 2]
logical_species <- logical_species[, c(-1, -2)]

# produce a data set with the number of ASVs in each (non-zero) OTU for each sample
# collect ASVs still present in the data after otu clustering, lulu curation , taxonomic assignment and firther filtering
sp_asv_otu <- subset(asv_otu, asv_otu$MOTU %in% rownames(species_MOTU))

all_asvs <- c()
for (row in 1:nrow(sp_asv_otu)){
  row_asvs <- grep("ASV", sp_asv_otu[row, 4:ncol(sp_asv_otu)], value = TRUE)
  row_asvs <- as.character(row_asvs)
  all_asvs <- c(all_asvs, row_asvs)
}

all_asvs <- unique(all_asvs)

# filter the imported asv table based on remaining asvs
sp_asv_presence_absence <- asv_count[asv_count$X %in% all_asvs, ]
rownames(sp_asv_presence_absence) <- sp_asv_presence_absence[, 1]
sp_asv_presence_absence <- sp_asv_presence_absence[, -1]
sp_asv_presence_absence[sp_asv_presence_absence > 0] <- 1

# filter the imported asv table based on samples remaning
sp_asv_presence_absence <- sp_asv_presence_absence[, colnames(sp_asv_presence_absence) %in% sample_names(rarefied_phylo)]

# we now have two presence/absence tables, one for OTUs and one for ASVs
# asv_presence_absence and presence_absence_otus
# summarize them for the number of ASVs/OTUs in each sample
sp_asv_richness <- apply(sp_asv_presence_absence, 2, sum)
sp_otu_richness <- apply(logical_species, 2, sum)

################# make mean species/ASV plots #########################


# obtain all sub-categories of chosen category
meta_cats <- sort(unique(sample_data(species_subset)$RetrievedYear))
class(meta_cats) <- "character"

# create empty data base
sum_summary <- data.frame(
  Year = c(),
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
  cat_subset <- subset_samples(species_subset, RetrievedYear == current_cat)
  # extract sample names from samples present in the subset
  cat_samples <- sample_names(cat_subset)

  # subset the samples to the ASV/OTU sums
  OTU_sum <- sp_otu_richness[cat_samples]
  ASV_sum <- sp_asv_richness[cat_samples]

  # calculate mean OTUs per sample within the sub-category 
  OTU_mean <- mean(OTU_sum, na.rm = TRUE)
  OTU_sd <- sd(OTU_sum)

  # calculate the mean ASV number per sample within the sub-category
  ASV_mean <- mean(ASV_sum, na.rm = TRUE)
  ASV_sd <- sd(ASV_sum)

  year <- rep(current_cat, length(OTU_sum))
  names(year) <- names(OTU_sum)

  # collect ASV_sum, OTU_sum, ASV_mean and OTU_mean in a data base
    output1 <- data.frame(
      Year = year,
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
species_OTU_mean_plot <- ggplot(sum_summary, aes(x = Year, y = OTU_sum)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
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

ggsave(plot = species_OTU_mean_plot, file = file.path(scatter_dir, "species_mean_sd.png"),
       dpi = 400)

  # plot the point plot
species_ASV_mean_plot <- ggplot(sum_summary, aes(x = Year, y = ASV_sum)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5) +
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

  ggsave(plot = species_ASV_mean_plot, file = file.path(scatter_dir, "species_mean_sd_ASVs.png"),
       dpi = 400)



################# make NMDS (OTU) plots #########################

nmds_presence_absence <- t(logical_species)
nmds_presence_absence <- as.matrix(nmds_presence_absence)

nmds <- metaMDS(nmds_presence_absence, distance = "jaccard")
# stressplot(nmds)
nmds_stress <- nmds$stress

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
data_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
data_scores <- transform(merge(sample_data(species_subset)[, 3], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(species_subset)[, 5], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(species_subset)[, 8], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
data_scores <- transform(merge(sample_data(species_subset)[, 9], data_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
species_nmds_year_seq <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  annotate("text", x = 0.4, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_year_seq, file = file.path(nmds_dir, "species_NMDS_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
species_nmds_year_unit <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  annotate("text", x = 0.4, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_year_unit, file = file.path(nmds_dir, "species_NMDS_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
species_nmds_year_frac <- ggplot(data_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  annotate("text", x = 0.4, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_year_frac, file = file.path(nmds_dir, "species_NMDS_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



################# make NMDS (ASV) plots #########################

nmds_asv <- t(sp_asv_presence_absence)

nmds <- metaMDS(nmds_asv, distance = "jaccard")
# stressplot(nmds)
nmds_stress <- nmds$stress

# following example from https://jkzorz.github.io/2019/06/06/NMDS.html

#extract NMDS scores (x and y coordinates)
asv_scores <- as.data.frame(scores(nmds)$sites)

# Source - https://stackoverflow.com/a/17376106
# Posted by thelatemail, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-09, License - CC BY-SA 3.0
asv_scores <- transform(merge(sample_data(species_subset)[, 3], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(species_subset)[, 5], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(species_subset)[, 8], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)
asv_scores <- transform(merge(sample_data(species_subset)[, 9], asv_scores, by = 0, all = TRUE), row.names = Row.names, Row.names = NULL)

# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
species_nmds_asv_year_seq <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_asv_year_seq, file = file.path(nmds_dir, "species_NMDS_ASV_RetrievedYear_SeqMethod.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
species_nmds_asv_year_unit <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = ARMSUnit, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_asv_year_unit, file = file.path(nmds_dir, "species_NMDS_ASV_RetrievedYear_ARMSUnit.png"),
       dpi = 400, width  = 8.76, height = 6.8)


# plot nmds plot with RetrievedYear and ARMSUnit
species_nmds_asv_year_frac <- ggplot(asv_scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  annotate("text", x = 0.75, y = 0.65, label = paste("2d stress =", round(nmds_stress, 3))) +
  scale_color_manual(values = gradient_colors) +
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

ggsave(plot = species_nmds_asv_year_frac, file = file.path(nmds_dir, "species_NMDS_ASV_RetrievedYear_FractionGroup.png"),
       dpi = 400, width  = 8.76, height = 6.8)



stop("End of script")

# downlaoded alien species subset and red list from artdatabanken, subset phyloseq to these lists
alien_species <- intersect(tax_table(rarefied_phylo)[, 6], alien_list$Vetenskapligt.namn)
red_list_species <- intersect(tax_table(rarefied_phylo)[, 6], red_list$Vetenskapligt.namn)

alien_phylo <- subset_taxa(rarefied_phylo, species %in% alien_species)
red_list_phylo <- subset_taxa(rarefied_phylo, species %in% red_list_species)

# plot the same graphs but for ambi group very sensitive to disturbance
# first merge genus and species to get species names
traits_species <- paste(traits_table$Genus, traits_table$Species)
traits_common <- intersect(tax_table(rarefied_phylo)[, 6], traits_species)
traits_phylo <- subset_taxa(rarefied_phylo, species %in% traits_common)

######### identify potential mismatches using species list from Artdatabanken ##########

# remove species associated with brackish water, forest, agricultural lands, urban environments,
# fjäll, wetlands and beaches
skagerrak_table <- skagerrak_table[grepl("\\(M\\)", skagerrak_table$Landskapstyp) &
                                  !grepl("\\((B|S|J|U|F|V|L|H)\\)", skagerrak_table$Landskapstyp), ]

# extract species from the data not present in the data from Artdatabanken 
# at artfakta.se/sok/taxa Geografi>Provins>Bohuslän, Geografi>Provins>Skagerrak and Ekologi>Landskapstyp>Hav
unexpected_species <- setdiff(tax_table(species_phylo)[, 6], skagerrak_table$Vetenskapligt.namn)
expected_species <- intersect(tax_table(species_phylo)[, 6], skagerrak_table$Vetenskapligt.namn)

expected_phylo <- subset_taxa(species_phylo, species %in% expected_species)


