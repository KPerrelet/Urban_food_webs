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
plot1 <- ggplot(fw.df, aes(x = Urbanization, y = Predator.frac)) + 
  geom_point(color = "goldenrod1") + 
  geom_smooth(method = "lm", fill = "goldenrod1", color = "black") + 
  labs(x = "", y = "Proportion of\npredators") +
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_blank(), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0, "cm")) 

plot2 <- ggplot(fw.df, aes(x = Urbanization, y = Link.Density)) + 
  geom_point(color = "#682d94ff") + 
  geom_smooth(method = "lm", fill = "#682d94ff", color = "black") + 
  labs(x = "", y = "Link density") +
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_blank(), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0, "cm")) 

plot3 <- ggplot(fw.df, aes(x = Urbanization, y = Frac.Cluster.BG)) + 
  geom_point(color = "mediumseagreen") + 
  geom_smooth(method = "lm", fill = "mediumseagreen", color = "black") + 
  labs(x = "Urbanization", y = "Aquatic-terrestrial\nmodules relative sizes") +
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0, "cm")) 

plot4 <- ggplot(fw.df, aes(x = Urbanization, y = Mean.ShortestPath)) + 
  geom_point(color = "royalblue2") + 
  geom_smooth(method = "lm", fill = "royalblue2", color = "black") + 
  labs(x = "", y = "Mean shortest distance\nbetween nodes") +
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_blank(), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.margin = margin(0, 0, 0, 0, "cm")) 

plot5 <- ggplot(fw.df, aes(x = Urbanization, y = Modularity)) + 
  geom_point(color = "grey50") + 
  geom_smooth(method = "lm", fill = "grey50", color = "black") + 
  labs(x = "Urbanization", y = "Modularity") +
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        axis.text.x =  element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title = element_text(size = 16), 
        plot.margin = margin(0, 0, 0, 0, "cm")) 

ggarrange(plot2, plot4, plot1, plot3, plot5, 
          nrow = 3, ncol = 2, align = "hv")

# Select relevant columns for correlation analysis
Modularity.df <- fw.df[, which(colnames(fw.df) %in% c("Urbanization",
                                                        "Link.Density",
                                                        "Predator.frac",
                                                        "Frac.Cluster.BG", 
                                                        "Mean.ShortestPath",
                                                        "Modularity"))] %>%
  relocate(Urbanization, Predator.frac, Link.Density, Frac.Cluster.BG, Mean.ShortestPath, Modularity)

# Compute correlation matrix for selected metrics
corr_mat <- round(cor(Modularity.df),2)
# Rename correlation matrix rows and columns for better readability
colnames(corr_mat) <- c("Urbanization", 
                        "Aquatic-terrestrial\nmodules relative sizes", 
                        "Proportion of\npredators", 
                        "Link density", 
                        "Mean shortest distance\nbetween nodes",
                        "Modularity")
rownames(corr_mat) <- c("Urbanization", 
                        "Aquatic-terrestrial\nmodules relative sizes", 
                        "Proportion of\npredators", 
                        "Link density", 
                        "Mean shortest distance\nbetween nodes",
                        "Modularity")

# Set diagonal values to 0 to avoid self-correlation visualization
diag(corr_mat) <- 0

# Generate correlation matrix visualization
ggcorrplot(corr_mat, method = "circle", type = "upper", insig = "blank", show.diag = T) +
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
