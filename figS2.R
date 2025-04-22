################################################################################
#################################### fig S2 ####################################
################################################################################

library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)

# Function to generate density distribution plots
make_distribution_plot <- function(df, var, label = NULL) {
  # Convert specified variable and its randomized counterpart into long format
  df.long <- pivot_longer(df[, colnames(df) %in% c(var, paste0(var, ".rand"))], 
                          cols = c(var, paste0(var, ".rand")))
  
  # Compute density for observed data
  d <- data.frame(x = density(df[, var])$x, y = density(df[, var])$y)
  # Filter density values to focus on the central 90%
  d$y[abs(cumsum(d$y * mean(diff(d$x))) - 0.5) > 0.45] <- 0
  
  # Repeat the same steps for randomized data
  d2 <- data.frame(x = density(df[, paste0(var, ".rand")])$x, y = density(df[, paste0(var, ".rand")])$y)
  d2$y[abs(cumsum(d2$y * mean(diff(d2$x))) - 0.5) > 0.45] <- 0
  
  # Create the density plot
  plot <- ggplot(df.long, aes(x = value, fill = name)) +
    geom_area(aes(x = x, y = y), fill = "#b2df8a", alpha = 1,
              data = d, inherit.aes = FALSE) +
    geom_area(aes(x = x, y = y), fill = "#1f78b4", alpha = 1,
              data = d2, inherit.aes = FALSE) +
    geom_density(alpha = 0, linetype = "solid") + 
    geom_histogram(aes(y = stat(density)), binwidth = (max(df.long$value) - min(df.long$value))/100, alpha = 0.3, color = "black", position = "identity") + 
    scale_fill_manual(values = c("#b2df8a", "#1f78b4"), 
                      labels = c("Observed", "Randomised")) + 
    labs(x = ifelse(is.null(label), var, label), y = "Density") + 
    theme_pubr(base_size = 18) + 
    theme(strip.background = element_blank(),
          panel.spacing = unit(1, "lines"), 
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 16), 
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 16)) +
    guides(fill = guide_legend(title = "Model", override.aes = list(alpha = 1), nrow = 2))
  
  return(plot)
}

# Load data
fw.df <- read.csv("foodweb_metrics.csv")

# Generate individual distribution plots for various metrics
Modularity.plot <- make_distribution_plot(fw.df, "Modularity")
Nestedness.plot <- make_distribution_plot(fw.df, "Nestedness")
Niche.plot <- make_distribution_plot(fw.df, "Nic.over.Con.Jac", label = "Niche overlap")
Skew.Degree.plot <- make_distribution_plot(fw.df, "Skew.Degree", label = "Degree skewness")
Mean.TL.plot <- make_distribution_plot(fw.df, "Mean.TL", label = "Mean trophic level")
Omnivory.plot <- make_distribution_plot(fw.df, "Omnivory")
Coherence.plot <- make_distribution_plot(fw.df, "Coherence")
MeanShortestPath.plot <- make_distribution_plot(fw.df, "Mean.ShortestPath", label = "Mean shortest path")
Frac.Cluster.BG.plot <- make_distribution_plot(fw.df, "Frac.Cluster.BG", label = "Aquatic-terrestrial module relative sizes")

# Extract the legend from one of the plots
# This allows for a shared legend across multiple plots
g <- ggplotGrob(Modularity.plot + theme(legend.position = "bottom"))
legend.grob <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

plot_grid(
  Nestedness.plot + theme(legend.position = "none"), 
  Modularity.plot + theme(legend.position = "none"), 
  Niche.plot + theme(legend.position = "none"), 
  Skew.Degree.plot + theme(legend.position = "none"),
  Mean.TL.plot + theme(legend.position = "none"), 
  Omnivory.plot + theme(legend.position = "none"), 
  Coherence.plot + theme(legend.position = "none"), 
  MeanShortestPath.plot + theme(legend.position = "none"),
  Frac.Cluster.BG.plot + theme(legend.position = "none"),
  legend.grob, 
  ncol = 2)
