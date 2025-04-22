################################################################################
##################################### fig 3 #################################### 
################################################################################

library(dplyr)
library(bruceR)
library(ggplot2)
library(ggpubr)

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

# Calculate the fraction of predator species within total species richness
fw.df$Predator.frac <- fw.df$Predator / fw.df$Sp.Richness

# Create a visualization of the predator fraction vs different covariates
frac.plot <- ggplot(fw.df, aes(y = Predator.frac)) + 
  stat_smooth(geom="line", aes(x = Urbanization, color = "#636363"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="ribbon", aes(x = Urbanization, fill = "#636363"), method = "lm", alpha = 0.3, size = 1) + 
  stat_smooth(geom="line", aes(x = Quantity, color = "#a6d96a"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="ribbon", aes(x = Quantity, fill = "#a6d96a"), method = "lm", alpha = 0.4, size = 1) + 
  stat_smooth(geom="line", aes(x = Quality, color = "#1b7837"), method = "lm", linetype = "dashed", alpha = 0.7, size = 1.3) +
  stat_smooth(geom="ribbon", aes(x = Quality, fill = "#1b7837"), method = "lm", alpha = 0.1, size = 1) +
  stat_smooth(geom="line", aes(x = Connectivity, color = "#92c5de"), method = "lm", alpha = 1, size = 1.3) +
  stat_smooth(geom="ribbon", aes(x = Connectivity, fill = "#92c5de"), method = "lm", alpha = 0.4, size = 1) + 
  stat_smooth(geom="line", aes(x = Complementarity, color = "#2166ac"), method = "lm", linetype = "dashed", alpha = 0.7, size = 1.3) +
  stat_smooth(geom="ribbon", aes(x = Complementarity, fill = "#2166ac"), method = "lm", alpha = 0.1, size = 1) +
  scale_color_manual(name = "",
                     values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                     labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  scale_fill_manual(name = "",
                    values = c("#2166ac", "#1b7837", "#636363", "#92c5de", "#a6d96a"),
                    labels = c("Complementarity", "Quality", "Urbanization", "Connectivity", "Quantity")) +
  labs(x = "", y = "Proportion of predators") +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        legend.position = "right",
        axis.text.x =  element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 18)) +
  guides(fill = guide_legend(nrow = 2), 
         color = guide_legend(nrow = 2))


# Standardize network structure and composition variables
variables.to.scale <- c("Sp.Richness", "Connectance", "Nestedness", "Modularity", "Nic.over.Con.Jac", "Skew.Degree", "Mean.TL", "Mean.Gen", "Omnivory", "Coherence")
fw.df[variables.to.scale] <- lapply(fw.df[variables.to.scale], scale)

# Perform linear regressions and store results
res <- data.frame()
for (response in variables.to.scale) {
  response.lm <- lm(fw.df[, response] ~ Predator.frac, data = fw.df)
  if (response %in% c("Sp.Richness", "Connectance", "Nestedness", "Modularity", "Nic.over.Con.Jac")) { # We create two groups (compositional and structural properties) for visualization purposes
    res <- rbind(res, data.frame(response = response, 
                                 coef = response.lm$coefficients[-1],
                                 # stderr = sqrt(diag(vcov(response.lm)))[-1], 
                                 CI2.5 = confint(response.lm)[2, 1], 
                                 CI97.5 = confint(response.lm)[2, 2], 
                                 pvalue = summary(response.lm)$coefficients[2, 4], 
                                 group = "Structure")) 
  }
  if (response %in% c("Skew.Degree", "Mean.TL", "Mean.Gen", "Omnivory", "Coherence")) {
    res <- rbind(res, data.frame(response = response, 
                                 coef = response.lm$coefficients[-1],
                                 # stderr = sqrt(diag(vcov(response.lm)))[-1], 
                                 CI2.5 = confint(response.lm)[2, 1], 
                                 CI97.5 = confint(response.lm)[2, 2], 
                                 pvalue = summary(response.lm)$coefficients[2, 4], 
                                 group = "Composition"))
  }
}

res$response <- paste0(res$response, "&", res$group)

# Create a coefficient plot
coef.plot <- ggplot(res, aes(x = coef, y = response)) +
  geom_tile(data = res[res$response %in% c(
    "Nestedness&Structure",
    "Connectance&Structure",
    "Coherence&Composition", 
    "Mean.Gen&Composition",
    "Skew.Degree&Composition"), ], 
    aes(y = response, x = 0, height = 1, width = Inf), fill =  "grey95", color = NA) +
  geom_vline(xintercept = 0) +
  geom_errorbar(aes(xmin = CI2.5, xmax = CI97.5), width = 0, linewidth = 1.2, position = position_dodge(width=0.7)) +
  geom_point(size = 3, position = position_dodge(width = 0.7)) +
  scale_y_discrete(
    limits = c("Nic.over.Con.Jac&Structure", 
               "Nestedness&Structure",
               "Modularity&Structure",
               "Connectance&Structure",
               "Sp.Richness&Structure",
               "Coherence&Composition", 
               "Omnivory&Composition",
               "Mean.Gen&Composition",
               "Mean.TL&Composition",
               "Skew.Degree&Composition"),
    labels = c("Niche overlap&Structure",
               "Nestedness&Structure",
               "Modularity&Structure",
               "Connectance&Structure",
               "Family richness&Structure",
               "Incoherence&Composition", 
               "Omnivory&Composition",
               " Mean generality&Composition",
               "Mean trophic\nlevel&Composition",
               "Node degree\nskewness&Composition")) +
  guides(y = ggh4x::guide_axis_nested(delim = "&")) + 
  labs(x = "Standardized coefficients and 95% CI", y = "") + 
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), 
        ggh4x.axis.nesttext.y = element_text(angle = 90, hjust = .5, vjust = 1))  

ggarrange(coef.plot, frac.plot, nrow = 1, common.legend = F,
          labels = c("A", "B"), 
          widths = c(1, 0.8),
          font.label = list(size = 22),
          legend = "bottom")
