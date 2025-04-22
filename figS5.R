################################################################################
#################################### fig S5 ####################################
################################################################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(bruceR)

# Load data
fw.df <- read.csv("foodweb_metrics.csv")
covariates.df <- read.csv("covariates.csv")
fw.df <- merge(fw.df, covariates.df, by = "site_no")

# Standardize and scale (0-1) the "Green50" variable to represent 'Quantity'
# Green50 is the number of "green" pixels in a 50 m buffer
# It gets larger if 1) their is more vegetation, and/or 2) if the pond gets larger
fw.df$Quantity <- scaler(fw.df$green50)

# Define a vector of selected quality-related covariates
covariate.quality <- c(
  "natural_banks",
  "marcophyte",
  "vegetation_height_sd_50",
  "vegetation_density_50"
)
PCA.quality.df <- fw.df[, colnames(fw.df) %in% covariate.quality]
# Perform Principal Component Analysis (PCA) on quality-related covariates
PCA.quality = prcomp(PCA.quality.df, center = T, scale. = T)
# Store the first principal component as the 'Quality' metric
fw.df$Quality <- scaler(predict(PCA.quality)[, 1])

# Aquatic connectivity is calculated as the mean distance to ponds in a 500 m buffer, so the larger, the less connected. 
# Thus, we take the opposite values, so that it aligns with terrestrial connectivitiy (in which, larger values = more connected landscapes)
fw.df$aquatic_connectivity <- fw.df$aquatic_connectivity * -1 
# Define a vector of connectivity-related covariates and repeat the PCA
covariate.connectivity <- c(
  "terrestrial_connectivity",
  "aquatic_connectivity",
  "distance_pond",
  "distance_forest"
)
PCA.connectivity.df <- fw.df[, colnames(fw.df) %in% covariate.connectivity]
PCA.connectivity.df = prcomp(PCA.connectivity.df, center = T, scale. = T)
fw.df$Connectivity <- scaler(predict(PCA.connectivity.df)[, 1]*-1)

# Complementarity is calculated using the joint entropy function of the landscapemetrics package
# high entropy = small diversity of values in co-occurence matrix
# low entropy = high diversity of values in co-occurence matrix
fw.df$Complementarity <- scaler(fw.df$complementarity)

# Define a vector of urban-related covariates and repeat the PCA
covariate.urban <- c(
  "overwarming",
  "no2",
  "population",
  "grey500_frac"
)
PCA.urban.df <- fw.df[, colnames(fw.df) %in% covariate.urban]
PCA.urban = prcomp(PCA.urban.df, center = T, scale. = T)
fw.df$Urbanization <- scaler(predict(PCA.urban)[, 1])

# Compute the proportion of predator species relative to total species richness
fw.df$Predator.frac <- fw.df$Predator / fw.df$Sp.Richness

plot1 <- ggplot(fw.df, aes(y = Basal)) + 
  stat_smooth(geom="ribbon", aes(x = Urbanization, fill = "#636363"), method = "lm", alpha = 0.1, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Quantity, fill = "#addd8e"), method = "lm", alpha = 0.1, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Quality, fill = "#31a354"), method = "lm", alpha = 0.1, size = 1) +
  stat_smooth(geom="ribbon", aes(x = Connectivity, fill = "#9ecae1"), method = "lm", alpha = 0.1, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Complementarity, fill = "#3182bd"), method = "lm", alpha = 0.1, size = 1) +
  stat_smooth(geom="line", aes(x = Urbanization, color = "#636363"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Quantity, color = "#addd8e"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Quality, color = "#31a354"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Connectivity, color = "#9ecae1"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Complementarity, color = "#3182bd"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  scale_color_manual(name = "",
                     values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                     labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  scale_fill_manual(name = "",
                    values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                    labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  labs(x = "", y = "Basal consumer richness") +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.position = "right",
        axis.text.x =  element_text(size = 12), 
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16)) 

plot2 <- ggplot(fw.df, aes(y = Predator)) + 
  stat_smooth(geom="ribbon", aes(x = Urbanization, fill = "#636363"), method = "lm", alpha = 0.3, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Quantity, fill = "#a6d96a"), method = "lm", alpha = 0.4, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Quality, fill = "#1b7837"), method = "lm", alpha = 0.1, size = 1) +
  stat_smooth(geom="ribbon", aes(x = Connectivity, fill = "#92c5de"), method = "lm", alpha = 0.4, size = 1) + 
  stat_smooth(geom="ribbon", aes(x = Complementarity, fill = "#2166ac"), method = "lm", alpha = 0.1, size = 1) +
  stat_smooth(geom="line", aes(x = Urbanization, color = "#636363"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Quantity, color = "#a6d96a"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Quality, color = "#1b7837"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Connectivity, color = "#92c5de"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="line", aes(x = Complementarity, color = "#2166ac"), method = "lm", linetype = "dashed", alpha = 1, size = 1.3) +
  scale_color_manual(name = "",
                     values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                     labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  scale_fill_manual(name = "",
                    values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                    labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  labs(x = "", y = "Predator richness") +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.position = "right",
        axis.text.x =  element_text(size = 12), 
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16)) 

ggarrange(plot1, plot2, ncol = 2, nrow = 1, common.legend = T)
