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
library(ggvenn)
library(vctrs)
library(eulerr)
library(glmmTMB)
library(ggpubr)
library(effects)


set.seed(1)  

perma_dir <- "results/figures/stats"
if(!dir.exists(perma_dir)) dir.create(perma_dir)

scatter_dir <- "results/figures/scatter_2019-2021"
if(!dir.exists(scatter_dir)) dir.create(scatter_dir)

filter_feed <- read.csv("output/final_taxonomy/TaxaExport_2026-04-26_11.59.56_filterfeeders.csv",
                        header = TRUE, sep = ";")

alien_list <- read.csv("output/final_taxonomy/artfakta-frammande.csv",
                        header = TRUE, sep = ";")
red_list <- read.csv("output/final_taxonomy/artfakta-re-cr-en-vu.csv",
                        header = TRUE, sep = ";")

# read phyloseq object
rarefied_phylo <- readRDS(file = file.path("output/final_taxonomy/rarefied_phyloseq.rds"))
rarefied_phylo <- subset_samples(rarefied_phylo, RetrievedYear %in% c("2019", "2020", "2021"))
sample_names <- rownames(sample_data(rarefied_phylo))

# read presence/absence tables and richness tables (summarized number of asvs/otus)
asv_presence_absence <- read.csv("output/stats/ASV_presence_absence.csv", 
                                sep = "\t", header = TRUE)
rownames(asv_presence_absence) <- asv_presence_absence$X
asv_presence_absence <- asv_presence_absence[colnames(asv_presence_absence) %in% sample_names]

otu_presence_absence <- read.csv("output/stats/OTU_presence_absence.csv", 
                                sep = "\t", header = TRUE)
rownames(otu_presence_absence) <- otu_presence_absence$X
otu_presence_absence <- otu_presence_absence[colnames(otu_presence_absence) %in% sample_names]

asv_richness <- read.csv("output/stats/ASV_richness.csv", 
                         sep = "\t", header = FALSE)
asv_richness <- asv_richness[asv_richness[, 1] %in% sample_names, ]
                
otu_richness <- read.csv("output/stats/OTU_richness.csv",  
                         sep = "\t", header = FALSE)
otu_richness <- otu_richness[otu_richness[, 1] %in% sample_names, ]


trans_otu <- as.data.frame(t(otu_presence_absence))
colnames(trans_otu) <- rownames(otu_presence_absence)

perma_df_otu <- merge(as.data.frame(sample_data(rarefied_phylo)), trans_otu, by = "row.names")
rownames(perma_df_otu) <- perma_df_otu[, 1]
perma_df_otu <- perma_df_otu[, -1]

trans_asv <- as.data.frame(t(asv_presence_absence))
colnames(trans_asv) <- rownames(asv_presence_absence)

perma_df_asv <- merge(as.data.frame(sample_data(rarefied_phylo)), trans_asv, by = "row.names")
rownames(perma_df_asv) <- perma_df_asv[, 1]
perma_df_asv <- perma_df_asv[, -1]

# factors to include
factors <- c("RetrievedYear", "FractionGroup", "ARMSUnit", "SeqMethod")

# factors to NOT include
not_factors <- c("MaterialSampleID", "Sequenced", "Fraction", "DeployedDate", "RetrievedDate")

# remove NOT factors and NOT samples
perma_df_asv <- perma_df_asv[!(perma_df_asv$ARMSUnit == "VH1" & perma_df_asv$FractionGroup == "MF100"), ]
perma_df_asv <- perma_df_asv[!(perma_df_asv$ARMSUnit == "VH2" & perma_df_asv$FractionGroup == "MF100"), ]
perma_df_asv <- perma_df_asv[!(perma_df_asv$ARMSUnit == "VH2" & perma_df_asv$FractionGroup == "MF500"), ]

# remove NOT factors and NOT samples
perma_df_otu <- perma_df_otu[!(perma_df_otu$ARMSUnit == "VH1" & perma_df_otu$FractionGroup == "MF100"), ]
perma_df_otu <- perma_df_otu[!(perma_df_otu$ARMSUnit == "VH2" & perma_df_otu$FractionGroup == "MF100"), ]
perma_df_otu <- perma_df_otu[!(perma_df_otu$ARMSUnit == "VH2" & perma_df_otu$FractionGroup == "MF500"), ]

UnitFraction <- as.character(interaction(perma_df_asv$FractionGroup, perma_df_asv$ARMSUnit))
perma_df_asv <- cbind(UnitFraction, perma_df_asv)

UnitFraction <- as.character(interaction(perma_df_otu$FractionGroup, perma_df_otu$ARMSUnit))
perma_df_otu <- cbind(UnitFraction, perma_df_otu)


# check if my clusters differ in dispersion or location. if differ in dispersion ---> not suitable for permanova
# check if my clusters differ in dispersion or location. if differ in dispersion ---> not suitable for permanova
perma_df_otu <- perma_df_otu[order(perma_df_otu$RetrievedYear), ]

otu_matrix <- perma_df_otu[, 11:ncol(perma_df_otu)]
otu_matrix <- as.matrix(otu_matrix)
storage.mode(otu_matrix) <- "numeric"
otu_nmds <- metaMDS(otu_matrix, distance = "jaccard")

otu_stress <- otu_nmds$stress

otu_scores <- as.data.frame(scores(otu_nmds)$sites)

otu_centroids <- data.frame(
    RetrievedYear = c("2019", "2020", "2021"),
    centroid_x = c(mean(otu_scores$NMDS1[perma_df_otu$RetrievedYear == "2019"]),
                   mean(otu_scores$NMDS1[perma_df_otu$RetrievedYear == "2020"]),
                   mean(otu_scores$NMDS1[perma_df_otu$RetrievedYear == "2021"])),
    centroid_y = c(mean(otu_scores$NMDS2[perma_df_otu$RetrievedYear == "2019"]),
                   mean(otu_scores$NMDS2[perma_df_otu$RetrievedYear == "2020"]),
                   mean(otu_scores$NMDS2[perma_df_otu$RetrievedYear == "2021"]))
)

n_2019 <- length(grep("2019", perma_df_otu$RetrievedYear))
n_2020 <- length(grep("2020", perma_df_otu$RetrievedYear))
n_2021 <- length(grep("2021", perma_df_otu$RetrievedYear))

plot_data <- data.frame(
    FractionGroup = perma_df_otu$FractionGroup,
    RetrievedYear = perma_df_otu$RetrievedYear,
    NMDS1 = (otu_scores$NMDS1),
    NMDS2 = (otu_scores$NMDS2),
    xend = c(rep(otu_centroids[1,2], n_2019), rep(otu_centroids[2,2], n_2020), rep(otu_centroids[3,2], n_2021)),
    yend = c(rep(otu_centroids[1,3], n_2019), rep(otu_centroids[2,3], n_2020), rep(otu_centroids[3,3], n_2021))
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
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  #stat_ellipse(geom = "polygon", alpha = 0.04, aes(group = RetrievedYear),
  #             color = "black", fill = "blue") +
  geom_point(data = otu_centroids, aes(x = centroid_x, y = centroid_y),
             color = "black", size = 2, shape = 7) + 
  geom_segment(data = plot_data, aes(x = NMDS1, y = NMDS2, 
               xend = xend, yend = yend), alpha = 0.5) +
  annotate("text", x = 1, y = 0.7, label = paste("2d stress =", round(otu_stress, 3))) +
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


ggsave(plot = distance_plot, file = file.path(perma_dir, "19_21_OTU_centroid_distance.png"),
       dpi = 400, width = 8, height = 8)

########### PERMANOVA ###########

# create distance matrix
otu_dist <- vegdist(otu_matrix, method = "jaccard")
 
# assumptions
otu_dispersion2 <- betadisper(otu_dist, group = perma_df_otu$RetrievedYear, type = "centroid")

plot(otu_dispersion2)

anova(otu_dispersion2)

ctrl <- how(within = Within(type = "free"),
    plots = Plots(strata = perma_df_otu$UnitFraction))

adonis2(otu_dist ~ RetrievedYear, 
        data = perma_df_otu,
        permutations = ctrl)


###################### ASV NMDS and PERMANOVA ###########################

# check if my clusters differ in dispersion or location. if differ in dispersion ---> not suitable for permanova
# check if my clusters differ in dispersion or location. if differ in dispersion ---> not suitable for permanova
perma_df_asv <- perma_df_asv[order(perma_df_asv$RetrievedYear), ]

asv_matrix <- perma_df_asv[, 11:ncol(perma_df_asv)]
asv_matrix <- as.matrix(asv_matrix)
storage.mode(asv_matrix) <- "numeric"
asv_nmds <- metaMDS(asv_matrix, distance = "jaccard")

asv_stress <- asv_nmds$stress

asv_scores <- as.data.frame(scores(asv_nmds)$sites)

asv_centroids <- data.frame(
    RetrievedYear = c("2019", "2020", "2021"),
    centroid_x = c(mean(asv_scores$NMDS1[perma_df_asv$RetrievedYear == "2019"]),
                   mean(asv_scores$NMDS1[perma_df_asv$RetrievedYear == "2020"]),
                   mean(asv_scores$NMDS1[perma_df_asv$RetrievedYear == "2021"])),
    centroid_y = c(mean(asv_scores$NMDS2[perma_df_asv$RetrievedYear == "2019"]),
                   mean(asv_scores$NMDS2[perma_df_asv$RetrievedYear == "2020"]),
                   mean(asv_scores$NMDS2[perma_df_asv$RetrievedYear == "2021"]))
)

n_2019 <- length(grep("2019", perma_df_asv$RetrievedYear))
n_2020 <- length(grep("2020", perma_df_asv$RetrievedYear))
n_2021 <- length(grep("2021", perma_df_asv$RetrievedYear))

plot_data <- data.frame(
    FractionGroup = perma_df_asv$FractionGroup,
    RetrievedYear = perma_df_asv$RetrievedYear,
    NMDS1 = (asv_scores$NMDS1),
    NMDS2 = (asv_scores$NMDS2),
    xend = c(rep(asv_centroids[1,2], n_2019), rep(asv_centroids[2,2], n_2020), rep(asv_centroids[3,2], n_2021)),
    yend = c(rep(asv_centroids[1,3], n_2019), rep(asv_centroids[2,3], n_2020), rep(asv_centroids[3,3], n_2021))
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
  geom_point(size = 3, aes(shape = FractionGroup, colour = RetrievedYear)) + 
  #stat_ellipse(geom = "polygon", alpha = 0.04, aes(group = RetrievedYear),
  #             color = "black", fill = "blue") +
  geom_point(data = asv_centroids, aes(x = centroid_x, y = centroid_y),
             color = "black", size = 2, shape = 7) + 
  geom_segment(data = plot_data, aes(x = NMDS1, y = NMDS2, 
               xend = xend, yend = yend), alpha = 0.5) +
  annotate("text", x = 0.75, y = 0.6, label = paste("2d stress =", round(asv_stress, 3))) +
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


ggsave(plot = distance_plot, file = file.path(perma_dir, "19_21_asv_centroid_distance.png"),
       dpi = 400, width = 8, height = 8)

########### PERMANOVA ###########

# create distance matrix
asv_dist <- vegdist(asv_matrix, method = "jaccard")
 
# assumptions
asv_dispersion2 <- betadisper(asv_dist, group = perma_df_asv$RetrievedYear, type = "centroid")

plot(asv_dispersion2)

anova(asv_dispersion2)

ctrl <- how(within = Within(type = "free"),
    plots = Plots(strata = perma_df_asv$UnitFraction))

adonis2(asv_dist ~ RetrievedYear, 
        data = perma_df_asv,
        permutations = ctrl)


###################### investigate OTU richness #########################

# otu_richness <- read.csv("output/stats/OTU_richness.csv", 
                         # sep = "\t", header = FALSE)

rownames(otu_richness) <- otu_richness[, 1]
otu_richness <- otu_richness[rownames(otu_richness) %in% rownames(perma_df_otu), ]

otu_richness <- merge(perma_df_otu[, c("UnitFraction", "SeqMethod", "RetrievedYear")], otu_richness, by = "row.names")

otu_richness <- otu_richness[order(otu_richness$UnitFraction, otu_richness$RetrievedYear), ]

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
gradient_colors <- colorRampPalette(gradient_colors)(6)


gradient_colors7 <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#E64BAF",
  "#F7923A",  # warm orange
  "#D6C229")  # clean golden yellow

gradient_colors7 <- colorRampPalette(gradient_colors7)(4)

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

ggsave(plot = glmm_lines, file = file.path(perma_dir, "2019_2021_OTU_linear_model.png"),
       dpi = 400)

# https://stats.stackexchange.com/questions/272633/how-to-decide-whether-to-set-reml-to-true-or-false/272654#272654
# It’s generally good to use REML, if it is available, when you are interested in the magnitude of the random effects variances, 
# but never when you are comparing models with different fixed effects via hypothesis tests or information-theoretic criteria such as AIC.
# to compute the likelihood, we must set REML = FALSE

# gaussian linear mixed model
glmm_gaussian <- lmer(V2 ~ RetrievedYear + (1 | UnitFraction), REML = FALSE, data = otu_richness)
summary(glmm_gaussian)

# are the residuals distributed normally?
qqnorm(resid(glmm_gaussian))
qqline(resid(glmm_gaussian))
shapiro.test(resid(glmm_gaussian))
# no, they are not, gausssian model is not suitable

# glmm with a poisson distribution
glmm_poisson <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = otu_richness, family = poisson)
summary(glmm_poisson)

# how are the residuals distributed? 
qqnorm(resid(glmm_poisson))
qqline(resid(glmm_poisson))

# glmm with a negative binomial distribution
glmm_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = otu_richness, family = nbinom2)
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
glmm3_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = otu_richness, family = nbinom2)
glmm3_sum <- summary(glmm3_nbiom)

# transform values back to non-log transformed. only valid for estimate and std error
glmm3_coeff <- glmm3_sum$coefficients[[1]]
exp(glmm3_coeff)


###################### investigate ASV richness #########################

# otu_richness <- read.csv("output/stats/OTU_richness.csv", 
                         # sep = "\t", header = FALSE)

rownames(asv_richness) <- asv_richness[, 1]
asv_richness <- asv_richness[rownames(asv_richness) %in% rownames(perma_df_asv), ]

asv_richness <- merge(perma_df_asv[, c("UnitFraction", "SeqMethod", "RetrievedYear")], asv_richness, by = "row.names")

asv_richness <- asv_richness[order(asv_richness$UnitFraction, asv_richness$RetrievedYear), ]

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
gradient_colors <- colorRampPalette(gradient_colors)(6)

asv_glmm_lines <- ggplot(asv_richness, aes(x = RetrievedYear, y = V2, group = UnitFraction, color = UnitFraction)) + 
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

ggsave(plot = asv_glmm_lines, file = file.path(perma_dir, "2019_2021_ASV_linear_model.png"),
       dpi = 400)

# https://stats.stackexchange.com/questions/272633/how-to-decide-whether-to-set-reml-to-true-or-false/272654#272654
# It’s generally good to use REML, if it is available, when you are interested in the magnitude of the random effects variances, 
# but never when you are comparing models with different fixed effects via hypothesis tests or information-theoretic criteria such as AIC.
# to compute the likelihood, we must set REML = FALSE

# gaussian linear mixed model
glmm_gaussian <- lmer(V2 ~ RetrievedYear + (1 | UnitFraction), REML = FALSE, data = asv_richness)
summary(glmm_gaussian)

# are the residuals distributed normally?
qqnorm(resid(glmm_gaussian))
qqline(resid(glmm_gaussian))
shapiro.test(resid(glmm_gaussian))
# no, they are not, gausssian model is not suitable

# glmm with a poisson distribution
glmm_poisson <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = asv_richness, family = poisson)
summary(glmm_poisson)

# how are the residuals distributed? 
qqnorm(resid(glmm_poisson))
qqline(resid(glmm_poisson))

# glmm with a negative binomial distribution
glmm_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = asv_richness, family = nbinom2)
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
glmm3_nbiom <- glmmTMB(V2 ~ RetrievedYear + (1 | UnitFraction), data = asv_richness, family = nbinom2)
glmm3_sum <- summary(glmm3_nbiom)

# transform values back to non-log transformed. only valid for estimate and std error
glmm3_coeff <- glmm3_sum$coefficients[[1]]
exp(glmm3_coeff)


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
        
      axis.text.x = element_text(hjust = 1, size = 10.5, angle = 45),
        
      axis.line.x = element_line(color = "#1e1e21", size = 0.5),
      axis.line.y = element_line(color = "#1e1e21", size = 0.5),
        
      axis.text = element_text(color = "#1e1e21"),
      axis.title = element_text(color = "#1e1e21", size = 12),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.ticks = element_line(color = "#1e1e21"),
      axis.ticks.length = unit(0.2, "cm"),
        
      panel.grid.major = element_line(color = "#babacc", size = 0.5, linetype = 3),
      panel.grid.minor = element_blank())

ggsave(plot = pres, file = "results/figures/2019_2021_species_presence_absence.png",
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

ggsave(plot = pres, file = "results/figures/2019_2021_genera_presence_absence.png",
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

ggsave(plot = pres, file = "results/figures/2019_2021_family_presence_absence.png",
       dpi = 400, height = 35, width = 7)


################ make venn diagram of species differences ###############

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A")  # warm orange)

# more chatgpt data transformation code
species_venn <- species_presence_absence %>%
  filter(Presence == 1)

species_venn <- species_venn %>%
  distinct(RetrievedYear, species)

species_venn <- species_venn %>%
  group_by(RetrievedYear) %>%
  summarise(species = list(as.character(species))) %>%
  deframe()

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A")  # clean golden yellow

fit <- euler(species_venn)

library(grid)
 
png(
  "euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


venn <- ggvenn(species_venn, 
fill_color = gradient_colors,
stroke_color = "#1e1e21",
stroke_size = 0.5) +
theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_blank(),
        
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
        
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
        
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

ggsave(plot = venn, file = "results/figures/2019_2021_venn.png",
       dpi = 400, height = 7, width = 7)


########## ASV and OTU euler plots, data transform with chatgpt ##########

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




otu_fit <- euler(otu_list)

library(grid)
 
png(
  "otu_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  otu_fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


asv_fit <- euler(asv_list)

library(grid)
 
png(
  "asv_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  asv_fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


####### investigate relationships between the years #######
all_years <- Reduce(intersect, species_venn)

unique_2019 <- setdiff(species_venn[["2019"]], unlist(species_venn[names(species_venn) != "2019"]))
unique_2020 <- setdiff(species_venn[["2020"]], unlist(species_venn[names(species_venn) != "2020"]))
unique_2021 <- setdiff(species_venn[["2021"]], unlist(species_venn[names(species_venn) != "2021"]))

phylo_2019 <- subset_taxa(species_phylo, species %in% unique_2019)
phylo_2020 <- subset_taxa(species_phylo, species %in% unique_2020)
phylo_2021 <- subset_taxa(species_phylo, species %in% unique_2021)

barplot_19 <- vec_count(as.vector(tax_table(phylo_2019)[, 1]))
barplot_19$year <- "2019"
barplot_20 <- vec_count(as.vector(tax_table(phylo_2020)[, 1]))
barplot_20$year <- "2020"
barplot_21 <- vec_count(as.vector(tax_table(phylo_2021)[, 1]))
barplot_21$year <- "2021"

df_all <- bind_rows(barplot_19, barplot_20, barplot_21)
colnames(df_all) <- c("phylum", "count", "year")

df_all <- as.data.frame(df_all)

df_all <- df_all %>%
  complete(phylum, year, fill = list(count = 0))

bar_sp <- ggplot(df_all, aes(x = phylum, y = count, fill = year)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() + 
    labs(x = "Phylum", y = "Count", fill = "Year") +
    scale_fill_manual(values = gradient_colors) +
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


######### model how the amount of filter feeders have changed #########

filter_phylo <- subset_taxa(species_phylo, species %in% filter_feed$Vetenskapligt.namn)
otu_presence_absence <- otu_presence_absence[rownames(otu_presence_absence) %in% rownames(otu_table(filter_phylo)), ]

otu_richness <- apply(otu_presence_absence, 2, sum)

# obtain all sub-categories of chosen category
meta_cats <- sort(unique(sample_data(filter_phylo)$RetrievedYear))
class(meta_cats) <- "character"

# create empty data base
sum_summary <- data.frame(
  Year = c(),
  OTU_sum = c()
)
mean_summary <- data.frame(
  OTU_mean = c(),
  OTU_sd = c()
)

for (cat in meta_cats){
  current_cat <- cat
  # create phyloseq object of the specific sub-category. e.g. 2019 or 2020
  cat_subset <- subset_samples(filter_phylo, RetrievedYear == current_cat)
  # extract sample names from samples present in the subset
  cat_samples <- sample_names(cat_subset)

  # subset the samples to the ASV/OTU sums
  OTU_sum <- otu_richness[cat_samples]

  # calculate mean OTUs per sample within the sub-category 
  OTU_mean <- mean(OTU_sum, na.rm = TRUE)
  OTU_sd <- sd(OTU_sum)

  year <- rep(current_cat, length(OTU_sum))
  names(year) <- names(OTU_sum)

  # collect ASV_sum, OTU_sum, ASV_mean and OTU_mean in a data base
    output1 <- data.frame(
      Year = year,
      OTU_sum = OTU_sum)
    output2 <- data.frame(
      OTU_mean = OTU_mean,
      OTU_sd = OTU_sd,
      row.names = as.character(cat))
  sum_summary <- rbind(sum_summary, output1)
  mean_summary <- rbind(mean_summary, output2)
}

mean_summary <- rownames_to_column(mean_summary, "Year") 


# investigate different plots (barplots, regular plots with trendlines, ggplot2)

# mean OTU plot with standard deviation as bars

# plot the point plot
OTU_mean_plot <- ggplot(sum_summary, aes(x = Year, y = OTU_sum)) + 
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

ggsave(plot = OTU_mean_plot, file = file.path(scatter_dir, "filter_feeders.png"),
       dpi = 400)



####### investigate richness ###########

# otu_richness <- read.csv("output/stats/OTU_richness.csv", 
                         # sep = "\t", header = FALSE)
filter_phylo <- subset_taxa(species_phylo, species %in% filter_feed$Vetenskapligt.namn)
otu_presence_absence <- otu_presence_absence[rownames(otu_presence_absence) %in% rownames(otu_table(filter_phylo)), ]

otu_richness <- apply(otu_presence_absence, 2, sum)

UnitFraction <- as.character(interaction(sample_data(filter_phylo)[[5]], sample_data(filter_phylo)[[9]]))

otu_richness <- merge(sample_data(filter_phylo)[, c(5, 8)], otu_richness, by = "row.names")
otu_richness <- cbind(UnitFraction, otu_richness)

otu_richness <- otu_richness[order(otu_richness$RetrievedYear), ]

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow
gradient_colors <- colorRampPalette(gradient_colors)(6)

glmm_lines <- ggplot(otu_richness, aes(x = RetrievedYear, y = y, group = UnitFraction, color = UnitFraction)) + 
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



ggsave(plot = glmm_lines, file = file.path(perma_dir, "2019_2021_filter_linear_model.png"),
       dpi = 400)

# https://stats.stackexchange.com/questions/272633/how-to-decide-whether-to-set-reml-to-true-or-false/272654#272654
# It’s generally good to use REML, if it is available, when you are interested in the magnitude of the random effects variances, 
# but never when you are comparing models with different fixed effects via hypothesis tests or information-theoretic criteria such as AIC.
# to compute the likelihood, we must set REML = FALSE


############### make euler plot of rare and alien species ############

gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A")

alien_species <- alien_list$Vetenskapligt.namn
red_species <- red_list$Vetenskapligt.namn

alien_venn <- lapply(species_venn, function(x) {
  x[x %in% alien_species]
})

red_venn <- lapply(species_venn, function(x) {
  x[x %in% red_species]
})


###### make alien venn/euler diagram ############
fit <- euler(alien_venn)

library(grid)
 
png(
  "alien_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


###### make red list venn/euler diagram ############
fit <- euler(red_venn)

library(grid)
 
png(
  "redlist_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()

############## do the same for all five years ###################
rarefied_phylo <- readRDS(file = file.path("output/final_taxonomy/rarefied_phyloseq.rds"))
sample_names <- rownames(sample_data(rarefied_phylo))

otu_presence_absence <- read.csv("output/stats/OTU_presence_absence.csv", 
                                sep = "\t", header = TRUE)
rownames(otu_presence_absence) <- otu_presence_absence$X
otu_presence_absence <- otu_presence_absence[colnames(otu_presence_absence) %in% sample_names]

trans_otu <- as.data.frame(t(otu_presence_absence))
colnames(trans_otu) <- rownames(otu_presence_absence)

perma_df_otu <- merge(as.data.frame(sample_data(rarefied_phylo)), trans_otu, by = "row.names")
rownames(perma_df_otu) <- perma_df_otu[, 1]
perma_df_otu <- perma_df_otu[, -1]

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



# more chatgpt data transformation code
species_venn <- species_presence_absence %>%
  filter(Presence == 1)

species_venn <- species_venn %>%
  distinct(RetrievedYear, species)

species_venn <- species_venn %>%
  group_by(RetrievedYear) %>%
  summarise(species = list(as.character(species))) %>%
  deframe()



gradient_colors <- c("#4C9FCD",  # soft blue (close to your original)
  "#6C4CCF",  # richer purple
  "#F7923A",  # warm orange
  "#E64BAF",  # vivid magenta
  "#D6C229")  # clean golden yellow

alien_venn <- lapply(species_venn, function(x) {
  x[x %in% alien_species]
})

red_venn <- lapply(species_venn, function(x) {
  x[x %in% red_species]
})


###### make alien venn/euler diagram ############
fit <- euler(alien_venn)

library(grid)
 
png(
  "all_alien_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


###### make red list venn/euler diagram ############
fit <- euler(red_venn)

library(grid)
 
png(
  "all_redlist_euler_plot.png",
  width = 7,
  height = 7,
  units = "in",
  res = 400,
  bg = "#fdfbf7",
)

grid.newpage()

# background
grid.rect(
  gp = gpar(fill = "#fdfbf7", col = NA)
)

# draw Euler plot
plot(
  fit,
  legend = list(side = "right", pch = 22, fontfamily = "sans", fontsize = "13", col = "#1e1e21", ncol  = 3),
  newpage = FALSE,
  quantities = list(
    type = c("counts", "percent"),
    fontsize = 13,
    fontfamily = "sans"
  ),
  fills = list(
    fill = gradient_colors,
    alpha = 0.6
  ),
  edges = list(
    col = "#1e1e21",
    lwd = 1
  ),
  labels = FALSE
)

dev.off()


all_red_venn <- ggvenn(red_venn, 
fill_color = gradient_colors,
stroke_color = "#1e1e21",
stroke_size = 0.5) +
theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_blank(),
        
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
        
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
        
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

ggsave(plot = all_red_venn, file = "results/figures/red_venn.png",
       dpi = 400, height = 7, width = 7)



all_alien_venn <- ggvenn(alien_venn, 
fill_color = gradient_colors,
stroke_color = "#1e1e21",
stroke_size = 0.5) +
theme_minimal() + 
  theme(plot.background = element_rect(fill = "#fdfbf7"),
      panel.background = element_blank(),
      panel.border = element_blank(),
        
      axis.text.x = element_blank(),
        
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
        
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
        
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

ggsave(plot = all_alien_venn, file = "results/figures/alien_venn.png",
       dpi = 400, height = 7, width = 7)

ggplot(data, aes(x=name, y=value)) + 
  geom_bar(stat = "identity")

