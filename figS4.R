################################################################################
##################################### fig 2 #################################### 
################################################################################

library(bruceR)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)

# Load data
fw.df <- read.csv("foodweb_metrics.csv")
covariates.df <- read.csv("covariates.csv")
fw.df <- merge(fw.df, covariates.df, by = "site_no")

# Define a vector of selected urban-related covariates
covariate.urban <- c(
  "overwarming",
  "no2",
  "population",
  "grey500_frac"
)
PCA.urban.df <- fw.df[, colnames(fw.df) %in% covariate.urban]
# Perform Principal Component Analysis (PCA) on urban-related covariates
PCA.urban = prcomp(PCA.urban.df, center = T, scale. = T)
# Store and scale (0-1) the first principal component as the 'Urban' metric
fw.df$Urbanization <- scaler(predict(PCA.urban)[, 1])

# Compute the proportion of predator species relative to total species richness
fw.df$Predator.frac <- fw.df$Predator / fw.df$Sp.Richness

# Create scatter plots of Urbanization vs several metrics (e.g., predator fraction, link density, relative size of the aquatic-terrestrial modules, etc.)
ggplot(fw.df, aes(x = Urbanization, y = N_Clust)) + 
  geom_point(color = "grey50") + 
  geom_smooth(method = "lm", fill = "grey50", color = "black") + 
  labs(x = "Urbanization", y = "Number of modules") +
  theme_pubr(base_size = 14) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        axis.text.x =  element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        plot.margin = margin(0, 0, 0, 0, "cm")) 
