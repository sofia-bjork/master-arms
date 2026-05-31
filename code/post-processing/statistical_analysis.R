#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(ggplot2)
library(tibble)
library(permute)
library(tidyverse)
library(lme4)
library(arm)
library(lmerTest)
library(betapart)
library(glmmTMB)
library(ggpubr)

set.seed(1)  

perma_dir <- "results/figures/stats"
if(!dir.exists(perma_dir)) dir.create(perma_dir)

# read phyloseq object
rarefied_phylo <- readRDS(file = file.path("output/final_taxonomy/rarefied_phyloseq.rds"))

# read presence/absence tables and richness tables (summarized number of asvs/otus)
asv_presence_absence <- read.csv("output/stats/ASV_presence_absence.csv", 
                                sep = "\t", header = TRUE)
otu_presence_absence <- read.csv("output/stats/OTU_presence_absence.csv", 
                                sep = "\t", header = TRUE)
asv_richness <- read.csv("output/stats/ASV_richness.csv", 
                         sep = "\t", header = FALSE)
otu_richness <- read.csv("output/stats/OTU_richness.csv", 
                         sep = "\t", header = FALSE)

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

# factors to include
factors <- c("RetrievedYear", "FractionGroup", "ARMSUnit", "SeqMethod")

# factors to NOT include
not_factors <- c("MaterialSampleID", "Sequenced", "Fraction", "DeployedDate", "RetrievedDate")

# samples to NOT inlcude (VH2, VH1 MF100)

# remove NOT factors and NOT samples
perma_df_otu <- perma_df_otu[, !(colnames(perma_df_otu) %in% not_factors)]
perma_df_otu <- perma_df_otu[perma_df_otu$ARMSUnit != "VH2", ]
perma_df_otu <- perma_df_otu[!(perma_df_otu$ARMSUnit == "VH1" & perma_df_otu$FractionGroup == "MF100"), ]

# remove NOT factors and NOT samples
perma_df_asv <- perma_df_asv[, !(colnames(perma_df_asv) %in% not_factors)]
perma_df_asv <- perma_df_asv[perma_df_asv$ARMSUnit != "VH2", ]
perma_df_asv <- perma_df_asv[!(perma_df_asv$ARMSUnit == "VH1" & perma_df_asv$FractionGroup == "MF100"), ]

UnitFraction <- as.character(interaction(perma_df_asv$FractionGroup, perma_df_asv$ARMSUnit))
perma_df_asv <- cbind(UnitFraction, perma_df_asv)
perma_df_asv <- perma_df_asv[, !(colnames(perma_df_asv) %in% c("ARMSUnit", "FractionGroup"))]

UnitFraction <- as.character(interaction(perma_df_otu$FractionGroup, perma_df_otu$ARMSUnit))
perma_df_otu <- cbind(UnitFraction, perma_df_otu)
perma_df_otu <- perma_df_otu[, !(colnames(perma_df_otu) %in% c("ARMSUnit", "FractionGroup"))]

# check if my clusters differ in dispersion or location. if differ in dispersion ---> not suitable for permanova
perma_df_otu <- perma_df_otu[order(perma_df_otu$SeqMethod), ]

otu_matrix <- perma_df_otu[, 4:ncol(perma_df_otu)]
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric"
otu_nmds <- metaMDS(otu_matrix, distance = "jaccard")

otu_stress <- otu_nmds$stress

otu_scores <- as.data.frame(scores(otu_nmds)$sites)

otu_centroids <- data.frame(
    SeqMethod = c("miseq", "novaseq"),
    centroid_x = c(mean(otu_scores$NMDS1[perma_df_otu$SeqMethod == "miseq"]),
                   mean(otu_scores$NMDS1[perma_df_otu$SeqMethod == "novaseq"])),
    centroid_y = c(mean(otu_scores$NMDS2[perma_df_otu$SeqMethod == "miseq"]),
                   mean(otu_scores$NMDS2[perma_df_otu$SeqMethod == "novaseq"]))
)

n_miseq <- length(grep("miseq", perma_df_otu$SeqMethod))
n_novaseq <- length(grep("novaseq", perma_df_otu$SeqMethod))

plot_data <- data.frame(
    SeqMethod = perma_df_otu$SeqMethod,
    RetrievedYear = perma_df_otu$RetrievedYear,
    NMDS1 = (otu_scores$NMDS1),
    NMDS2 = (otu_scores$NMDS2),
    xend = c(rep(otu_centroids[1,2], n_miseq), rep(otu_centroids[2,2], n_novaseq)),
    yend = c(rep(otu_centroids[1,3], n_miseq), rep(otu_centroids[2,3], n_novaseq))
)


# create gradient of colors for phylum graph
gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
# show_col(gradient_colors)

# plot nmds plot with RetrievedYear and SeqMethod
distance_plot <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 3, aes(shape = SeqMethod, colour = RetrievedYear)) + 
  #stat_ellipse(geom = "polygon", alpha = 0.04, aes(group = SeqMethod),
  #             color = "black", fill = "blue") +
  geom_point(data = otu_centroids, aes(x = centroid_x, y = centroid_y),
             color = "black", size = 2, shape = 7) + 
  geom_segment(data = plot_data, aes(x = NMDS1, y = NMDS2, 
               xend = xend, yend = yend), alpha = 0.5) +
  annotate("text", x = 0.4, y = 0.65, label = paste("2d stress =", round(otu_stress, 3))) +
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

ggsave(plot = distance_plot, file = file.path(perma_dir, "OTU_centroid_distance.png"),
       dpi = 400, width = 8, height = 8)

########### PERMANOVA ###########

# create distance matrix
otu_dist <- vegdist(otu_matrix, method = "jaccard")

# assumptions
otu_dispersion1 <- betadisper(otu_dist, group = perma_df_otu$SeqMethod, type = "centroid")
otu_dispersion2 <- betadisper(otu_dist, group = perma_df_otu$RetrievedYear, type = "centroid")

plot(otu_dispersion1)
plot(otu_dispersion2)

anova(otu_dispersion1)
anova(otu_dispersion2)

ctrl <- how(within = Within(type = "free"),
    plots = Plots(strata = perma_df_otu$UnitFraction))

adonis2(otu_dist ~ SeqMethod, 
        data = perma_df_otu,
        permutations = ctrl)

adonis2(otu_dist ~ RetrievedYear, 
        data = perma_df_otu,
        permutations = ctrl)


######### investigate nestedness and turnover #############
perma_df_otu <- perma_df_otu[order(perma_df_otu$RetrievedYear, perma_df_otu$UnitFraction), ]

beta_2019 <- perma_df_otu[perma_df_otu$RetrievedYear == "2019", -c(2, 3)]
rownames(beta_2019) <- beta_2019$UnitFraction
beta_2019 <- beta_2019[, -1]
beta_2019 <- as.matrix(beta_2019)
storage.mode(beta_2019) <- "numeric"
# beta_2019 <- apply(beta_2019, 2, sum)
# beta_2019[beta_2019 > 0] <- 1

beta_2020 <- perma_df_otu[perma_df_otu$RetrievedYear == "2020", -c(2, 3)]
rownames(beta_2020) <- beta_2020$UnitFraction
beta_2020 <- beta_2020[, -1]
beta_2020 <- as.matrix(beta_2020)
storage.mode(beta_2020) <- "numeric"

beta_2021 <- perma_df_otu[perma_df_otu$RetrievedYear == "2021", -c(2, 3)]
rownames(beta_2021) <- beta_2021$UnitFraction
beta_2021 <- beta_2021[, -1]
beta_2021 <- as.matrix(beta_2021)
storage.mode(beta_2021) <- "numeric"

beta_2022 <- perma_df_otu[perma_df_otu$RetrievedYear == "2022", -c(2, 3)]
rownames(beta_2022) <- beta_2022$UnitFraction
beta_2022 <- beta_2022[, -1]
beta_2022 <- as.matrix(beta_2022)
storage.mode(beta_2022) <- "numeric"

beta_2023 <- perma_df_otu[perma_df_otu$RetrievedYear == "2023", -c(2, 3)]
rownames(beta_2023) <- beta_2023$UnitFraction
beta_2023 <- beta_2023[, -1]
beta_2023 <- as.matrix(beta_2023)
storage.mode(beta_2023) <- "numeric"

beta.temp(beta_2022, beta_2023, index.family = "jaccard")

###### heatmap of OTUs present each year ###################

sum_2019 <- colSums(beta_2019)
sum_2019 <- ifelse(sum_2019 > 0, 1, 0)

sum_2020 <- colSums(beta_2020)
sum_2020 <- ifelse(sum_2020 > 0, 1, 0)

sum_2021 <- colSums(beta_2021)
sum_2021 <- ifelse(sum_2021 > 0, 1, 0)

sum_2022 <- colSums(beta_2022)
sum_2022 <- ifelse(sum_2022 > 0, 1, 0)

sum_2023 <- colSums(beta_2023)
sum_2023 <- ifelse(sum_2023 > 0, 1, 0)

heatmap_data <- data.frame(
    Sum2019 = sum_2019,
    Sum2020 = sum_2020,
    Sum2021 = sum_2021,
    Sum2022 = sum_2022,
    Sum2023 = sum_2023)
heatmap_data <- t(heatmap_data)

heatmap_data <- as.matrix(heatmap_data)
storage.mode(heatmap_data) <- "numeric"

heatmap(heatmap_data)
plot_heatmap(rarefied_phylo)



###################### investigate OTU richness #########################

otu_richness <- read.csv("output/stats/OTU_richness.csv", 
                         sep = "\t", header = FALSE)

rownames(otu_richness) <- otu_richness[, 1]
otu_richness <- otu_richness[rownames(otu_richness) %in% rownames(perma_df_otu), ]

otu_richness <- merge(perma_df_otu[, 1:3], otu_richness, by = "row.names")

otu_richness <- otu_richness[order(otu_richness$UnitFraction, otu_richness$RetrievedYear), ]

glmm_lines <- ggplot(otu_richness, aes(x = RetrievedYear, y = V2, group = UnitFraction, color = UnitFraction)) + 
  geom_point() +
  geom_line() +
  scale_color_manual(values = gradient_colors) +
  theme_minimal() + 
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

ggsave(plot = glmm_lines, file = file.path(perma_dir, "OTU_linear_model.png"),
       dpi = 400)

# gaussian linear mixed model
glmm_gaussian <- lmer(V2 ~ SeqMethod + (1 | UnitFraction), REML = FALSE, data = otu_richness)
summary(glmm_gaussian)

# are the residuals distributed normally?
qqnorm(resid(glmm_gaussian))
qqline(resid(glmm_gaussian))
shapiro.test(resid(glmm_gaussian))
# no, they are not, gausssian model is not suitable

# glmm with a poisson distribution
glmm_poisson <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = otu_richness, family = poisson)
summary(glmm_poisson)

# how are the residuals distributed? 
qqnorm(resid(glmm_poisson))
qqline(resid(glmm_poisson))

# glmm with a negative binomial distribution
glmm_nbiom <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = otu_richness, family = nbinom2)
summary(glmm_nbiom)

# how are the residuals distributed?
glmm1ModelOutputs<-data.frame(Fitted=fitted(glmm_nbiom),
                  Residuals=resid(glmm_nbiom))

p5<-ggplot(glmm1ModelOutputs)+
    geom_point(aes(x=Fitted,y=Residuals))+
    theme_classic()+
    labs(y="Residuals",x="Fitted Values")

p6<-ggplot(glmm1ModelOutputs) +
    stat_qq(aes(sample=Residuals))+
    stat_qq_line(aes(sample=Residuals))+
    theme_classic()+
    labs(y="Sample Quartiles",x="Theoretical Quartiles")

# compare methods using AIC
AIC(glmm_gaussian, glmm_poisson)
AIC(glmm_gaussian, glmm_nbiom)
AIC(glmm_poisson, glmm_nbiom)

# glmm with a negative binomial distribution
glmm2_nbiom <- glmmTMB(V2 ~ SeqMethod + RetrievedYear + (1 | UnitFraction), data = otu_richness, family = nbinom2)
glmm2_sum <- summary(glmm2_nbiom)

glmm3_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = otu_richness, family = nbinom2)
glmm3_sum <- summary(glmm3_nbiom)

glmm4_nbiom <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = otu_richness, family = nbinom2)
glmm4_sum <- summary(glmm4_nbiom)

AIC(glmm2_nbiom, glmm3_nbiom, glmm4_nbiom)

# transform values back to non-log transformed. only valid for estimate and std error
glmm2_coeff <- glmm2_sum$coefficients[[1]]
exp(glmm2_coeff)
glmm3_coeff <- glmm3_sum$coefficients[[1]]
exp(glmm3_coeff)
glmm4_coeff <- glmm4_sum$coefficients[[1]]
exp(glmm4_coeff)

###################### make glmm of the years 2021, 2022 and 2023 ###############

subset_richness <- otu_richness[otu_richness$RetrievedYear %in% c("2021", "2022", "2023"), ]

glmm5_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = subset_richness, family = nbinom2)
glmm5_sum <- summary(glmm5_nbiom)

glmm6_nbiom <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = subset_richness, family = nbinom2)
glmm6_sum <- summary(glmm6_nbiom)

###################### investigate batch effects in OTU richness #########################
subset_richness <- otu_richness[otu_richness$RetrievedYear %in% c("2019", "2020", "2021"), ]

rownames(subset_richness) <- subset_richness$Row.names

batch_otu_richness <- merge(as.data.frame(sample_data(rarefied_phylo)[, 2]), subset_richness, by = "row.names")
batch_otu_richness[order(batch_otu_richness$Sequenced, decreasing = TRUE), ]


batch1_nbiom <- glmmTMB(V2 ~ RetrievedYear + Sequenced + (1 | UnitFraction), data = batch_otu_richness, family = nbinom2)
batch1_sum <- summary(batch1_nbiom)

batch2_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = batch_otu_richness, family = nbinom2)
batch2_sum <- summary(batch2_nbiom)

batch3_nbiom <- glmmTMB(V2 ~ Sequenced + (1 | UnitFraction), data = batch_otu_richness, family = nbinom2)
batch3_sum <- summary(batch3_nbiom)


###################### investigate ASV richness #########################

asv_richness <- read.csv("output/stats/ASV_richness.csv", 
                         sep = "\t", header = FALSE)

rownames(asv_richness) <- asv_richness[, 1]
asv_richness <- asv_richness[rownames(asv_richness) %in% rownames(asv_richness), ]

asv_richness <- merge(perma_df_asv[, 1:3], asv_richness, by = "row.names")

asv_richness <- asv_richness[order(asv_richness$UnitFraction, asv_richness$RetrievedYear), ]

glmm_lines <- ggplot(asv_richness, aes(x = RetrievedYear, y = V2, group = UnitFraction, color = UnitFraction)) + 
  geom_point() +
  geom_line() +
  scale_color_manual(values = gradient_colors) +
  theme_minimal() + 
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

ggsave(plot = glmm_lines, file = file.path(perma_dir, "ASV_linear_model.png"),
       dpi = 400)

# gaussian linear mixed model
glmm_gaussian <- lmer(V2 ~ SeqMethod + (1 | UnitFraction), REML = FALSE, data = asv_richness)
summary(glmm_gaussian)

# are the residuals distributed normally?
qqnorm(resid(glmm_gaussian))
qqline(resid(glmm_gaussian))
shapiro.test(resid(glmm_gaussian))

# glmm with a poisson distribution
glmm_poisson <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = asv_richness, family = poisson)
summary(glmm_poisson)

# how are the residuals distributed? 
qqnorm(resid(glmm_poisson))
qqline(resid(glmm_poisson))

# glmm with a negative binomial distribution
glmm_nbiom <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = asv_richness, family = nbinom2)
summary(glmm_nbiom)

# how are the residuals distributed?
glmm1ModelOutputs<-data.frame(Fitted=fitted(glmm_nbiom),
                  Residuals=resid(glmm_nbiom))

p5<-ggplot(glmm1ModelOutputs)+
    geom_point(aes(x=Fitted,y=Residuals))+
    theme_classic()+
    labs(y="Residuals",x="Fitted Values")

p6<-ggplot(glmm1ModelOutputs) +
    stat_qq(aes(sample=Residuals))+
    stat_qq_line(aes(sample=Residuals))+
    theme_classic()+
    labs(y="Sample Quartiles",x="Theoretical Quartiles")

# compare methods using AIC
AIC(glmm_gaussian, glmm_poisson)
AIC(glmm_gaussian, glmm_nbiom)
AIC(glmm_poisson, glmm_nbiom)

# glmm with a negative binomial distribution
glmm2_nbiom <- glmmTMB(V2 ~ SeqMethod + RetrievedYear + (1 | UnitFraction), data = asv_richness, family = nbinom2)
glmm2_sum <- summary(glmm2_nbiom)

glmm3_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = asv_richness, family = nbinom2)
glmm3_sum <- summary(glmm3_nbiom)

glmm4_nbiom <- glmmTMB(V2 ~ SeqMethod + (1 | UnitFraction), data = asv_richness, family = nbinom2)
glmm4_sum <- summary(glmm4_nbiom)

# transform values back to non-log transformed. only valid for estimate and std error
glmm2_coeff <- glmm2_sum$coefficients[[1]]
exp(glmm2_coeff)
glmm3_coeff <- glmm3_sum$coefficients[[1]]
exp(glmm3_coeff)
glmm4_coeff <- glmm4_sum$coefficients[[1]]
exp(glmm4_coeff)

###################### make glmm of the years 2021, 2022 and 2023 ###############

asv_subset_richness <- asv_richness[asv_richness$RetrievedYear %in% c("2021", "2022", "2023"), ]
glmm5_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = asv_subset_richness, family = nbinom2)
glmm5_sum <- summary(glmm5_nbiom)

########### plot species presence/absence 2019-2021 #############

species_phylo <- subset_taxa(rarefied_phylo, !(is.na(species)))
species_phylo <- prune_samples(sample_names(species_phylo) %in% rownames(perma_df_otu), species_phylo)
species_phylo <- prune_taxa(taxa_sums(species_phylo) > 0, species_phylo)
species_melt <- psmelt(species_phylo)

species_presence_absence <- species_melt %>%
  filter(!is.na(species)) %>%
  mutate(Presence = ifelse(Abundance > 0, 1, 0)) %>%
  group_by(RetrievedYear, species) %>%
  summarise(Presence = max(Presence), .groups = "drop")

species_presence_absence <- species_presence_absence %>%
  group_by(species) %>%
  mutate(freq = sum(Presence)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, freq))


pres <- ggplot(species_presence_absence, aes(x = factor(RetrievedYear), y = species, fill = factor(Presence))) + 
  geom_tile(colour = "#fdfbf7") +
  scale_fill_manual(values = c("0" = "#fdfbf7", "1" = "#1e1e21")) +
  labs(x = "Year", y = "Species", fill = "") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_text(hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank())

ggsave(plot = pres, file = "results/figures/species_presence_absence.png",
       dpi = 400, height = 35, width = 7)


########### plot genera presence/absence 2019-2021 #############

genera_phylo <- subset_taxa(rarefied_phylo, !(is.na(genus)))
genera_phylo <- prune_samples(sample_names(genera_phylo) %in% rownames(perma_df_otu), genera_phylo)
genera_phylo <- prune_taxa(taxa_sums(genera_phylo) > 0, genera_phylo)
genera_melt <- psmelt(genera_phylo)

genera_presence_absence <- genera_melt %>%
  filter(!is.na(genus)) %>%
  mutate(Presence = ifelse(Abundance > 0, 1, 0)) %>%
  group_by(RetrievedYear, genus) %>%
  summarise(Presence = max(Presence), .groups = "drop")

genera_presence_absence <- genera_presence_absence %>%
  group_by(genus) %>%
  mutate(freq = sum(Presence)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(genus, freq))


pres <- ggplot(genera_presence_absence, aes(x = factor(RetrievedYear), y = genus, fill = factor(Presence))) + 
  geom_tile(colour = "#fdfbf7") +
  scale_fill_manual(values = c("0" = "#fdfbf7", "1" = "#1e1e21")) +
  labs(x = "Year", y = "Genus", fill = "") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_text(hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank())

ggsave(plot = pres, file = "results/figures/genera_presence_absence.png",
       dpi = 400, height = 35, width = 7)


########### plot family presence/absence 2019-2021 #############

family_phylo <- subset_taxa(rarefied_phylo, !(is.na(family)))
family_phylo <- prune_samples(sample_names(family_phylo) %in% rownames(perma_df_otu), family_phylo)
family_phylo <- prune_taxa(taxa_sums(family_phylo) > 0, family_phylo)
family_melt <- psmelt(family_phylo)

family_presence_absence <- family_melt %>%
  filter(!is.na(family)) %>%
  mutate(Presence = ifelse(Abundance > 0, 1, 0)) %>%
  group_by(RetrievedYear, family) %>%
  summarise(Presence = max(Presence), .groups = "drop")

family_presence_absence <- family_presence_absence %>%
  group_by(family) %>%
  mutate(freq = sum(Presence)) %>%
  ungroup() %>%
  mutate(species = fct_reorder(family, freq))


pres <- ggplot(family_presence_absence, aes(x = factor(RetrievedYear), y = family, fill = factor(Presence))) + 
  geom_tile(colour = "#fdfbf7") +
  scale_fill_manual(values = c("0" = "#fdfbf7", "1" = "#1e1e21")) +
  labs(x = "Year", y = "Family", fill = "") +
  theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_text(hjust = 1, size = 10.5),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank())

ggsave(plot = pres, file = "results/figures/family_presence_absence.png",
       dpi = 400, height = 35, width = 7)


############ make lists of asvs and otus present year by year ################
# (data formatting done with help from chat gpt)

otu_list <- perma_df_otu %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = starts_with("MOTU"),
    names_to = "MOTU",
    values_to = "Present"
  ) %>%
  filter(Present == 1) %>%
  distinct(RetrievedYear, MOTU) %>%
  group_by(RetrievedYear) %>%
  summarise(OTU = list(sort(MOTU)), .groups = "drop") %>%
  deframe()

asv_list <- perma_df_asv %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = starts_with("ASV"),
    names_to = "ASV",
    values_to = "Present"
  ) %>%
  filter(Present == 1) %>%
  distinct(RetrievedYear, ASV) %>%
  group_by(RetrievedYear) %>%
  summarise(ASV = list(sort(ASV)), .groups = "drop") %>%
  deframe()

  # more chatgpt data transformation code
species_venn <- species_presence_absence %>%
  filter(Presence == 1)

species_venn <- species_venn %>%
  distinct(RetrievedYear, species)

species_venn <- species_venn %>%
  group_by(RetrievedYear) %>%
  summarise(species = list(as.character(species))) %>%
  deframe()

# number of shared asvs, otus and species

shared_asvs <- Reduce(intersect, asv_list)
shared_otus <- Reduce(intersect, otu_list)
shared_species <- Reduce(intersect, species_venn)

length(shared_asvs)
length(shared_otus)
length(shared_species)

# percentage of reads from shared species
sum(sample_sums(prune_taxa(shared_otus, rarefied_phylo))) / sum(sample_sums(rarefied_phylo))

# percentage of reads from shared species
sum(sample_sums(subset_taxa(rarefied_phylo, species %in% shared_species))) / sum(sample_sums(rarefied_phylo))

# number of unique asvs and otus each year, across all present samples
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

otu_list <- perma_df_otu %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = starts_with("MOTU"),
    names_to = "MOTU",
    values_to = "Present"
  ) %>%
  filter(Present == 1) %>%
  distinct(RetrievedYear, MOTU) %>%
  group_by(RetrievedYear) %>%
  summarise(OTU = list(sort(MOTU)), .groups = "drop") %>%
  deframe()

asv_list <- perma_df_asv %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = starts_with("ASV"),
    names_to = "ASV",
    values_to = "Present"
  ) %>%
  filter(Present == 1) %>%
  distinct(RetrievedYear, ASV) %>%
  group_by(RetrievedYear) %>%
  summarise(ASV = list(sort(ASV)), .groups = "drop") %>%
  deframe()


lengths(asv_list)
lengths(otu_list)
