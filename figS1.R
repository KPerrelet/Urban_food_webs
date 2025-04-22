################################################################################
#################################### fig S1 ####################################
################################################################################

library(ggbiplot)
library(dplyr)
library(ggcorrplot)
library(ggpubr)

# Load data
fw.df <- read.csv("foodweb_metrics.csv")
covariates.df <- read.csv("covariates.csv")
fw.df <- merge(fw.df, covariates.df, by = "site_no")

# Standardize the "Green50" variable to represent 'Quantity'
# Green50 is the number of "green" pixels in a 50 m buffer
# It gets larger if 1) their is more vegetation, and/or 2) if the pond gets larger
fw.df$Quantity <- scale(fw.df$green50)

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
# Create a biplot for quality
plot1 <- ggbiplot(PCA.quality,
                  point.size	= 2,
                  alpha = 0.5,
                  varname.size = 5) +
  labs(title = "Quality", x = "PC1", y = "PC2") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=14)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
# Store the first principal component as the 'Quality' metric
fw.df$Quality <- scale(predict(PCA.quality)[, 1])


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
plot2 <- ggbiplot(PCA.connectivity.df,
                  point.size	= 2,
                  alpha = 0.5,
                  varname.size = 5) +
  labs(title = "Connectivity", x = "PC1", y = "PC2") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=14)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
fw.df$Connectivity <- scale(predict(PCA.connectivity.df)[, 1]*-1)

# Complementarity is calculated using the joint entropy function of the landscapemetrics package
# high entropy = small diversity of values in co-occurence matrix
# low entropy = high diversity of values in co-occurence matrix
fw.df$Complementarity <- scale(fw.df$complementarity)

# Define a vector of urban-related covariates and repeat the PCA
covariate.urban <- c(
  "overwarming",
  "no2",
  "population",
  "grey500_frac"
)
PCA.urban.df <- fw.df[, colnames(fw.df) %in% covariate.urban]
PCA.urban = prcomp(PCA.urban.df, center = T, scale. = T)
plot3 <- ggbiplot(PCA.urban,
                  point.size	= 2,
                  alpha = 0.5,
                  varname.size = 5) +
  labs(title = "Urbanization", x = "PC1", y = "PC2") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size=14)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
fw.df$Urbanization <- scale(predict(PCA.urban)[, 1])

covariates.df <- fw.df[, which(colnames(fw.df) %in% c("Urbanization",
                                                        "Quality",
                                                        "Quantity",
                                                        "Connectivity",
                                                        "Complementarity"))] %>%
  relocate(Urbanization, Quantity, Quality, Connectivity, Complementarity)

corr_mat <- round(cor(covariates.df),2)
p.mat <- cor_pmat(covariates.df)
plot4 <- ggcorrplot(corr_mat, method = "circle", type = "upper", insig = "blank", show.diag = T, p.mat = p.mat) +
  labs(x = "", y = "") + 
  theme_void() + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text.y = element_text(size = 18, hjust = 1),
        axis.title = element_text(size = 18), 
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_fill_gradient2(low = "red", 
                       high = "blue", 
                       name = "Correlation")

ggarrange(plot1, plot2, plot3, plot4, 
          nrow = 2, ncol = 2,
          # align = "hv", 
          labels = c("A", "B", "C", "D"))
