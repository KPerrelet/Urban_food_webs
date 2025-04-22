################################################################################
#################################### fig S3 ####################################
################################################################################

library(piecewiseSEM)
library(igraph)
library(qgraph)

# Function to plot SEM results
plot_sem <- function(sem, res_labels, pred_labels, size = c(43, 30), edge.cex = 2.2) {
  # Extract coefficients from SEM summary
  sem.df <- summary(sem, .progressBar = F)$coefficients
  
  # Define edge styles based on p-values
  sem.df$style <- ifelse(sem.df$P.Value < 0.1, ifelse(sem.df$P.Value < 0.05, "solid", "dashed"), "solid")
  sem.df$color <- ifelse(sem.df$P.Value < 0.1, ifelse(sem.df$P.Value < 0.05, "black", "black"), "grey")
  
  # Assign edge labels only to significant paths
  sem.df[sem.df$color == "grey", ]$Std.Estimate <- NA
  sem.df[sem.df$color == "black", ]$Std.Estimate <- round(sem.df[sem.df$color == "black", ]$Std.Estimate, 3)
  
  # Assign colors based on effect direction (red for negative, blue for positive)
  sem.df[sem.df$color == "black", ]$color <- ifelse(sem.df[sem.df$color == "black", ]$Std.Estimate < 0, "#d7301f", "#0570b0")
  
  # Background colors for edges
  sem.df$color_bg <- NA
  sem.df[sem.df$color == "grey", ]$color_bg <- NA
  sem.df[sem.df$color == "#d7301f", ]$color_bg <- "#fddbc7"
  sem.df[sem.df$color == "#0570b0", ]$color_bg <- "#d1e5f0"
  
  # Create a directed acyclic graph (DAG) from SEM
  g <- graph_from_adjacency_matrix(getDAG(sem))
  EL <- as_edgelist(g)
  
  # Order data frame according to edge list
  sem.df <- sem.df[order(match(sem.df$Response, EL[, 2])), ]
  sem.df <- sem.df[order(match(sem.df$Predictor, EL[, 1])), ]
  
  # Define layout coordinates for node positioning
  if (length(res_labels) == 1) {
    if (res_labels %in% c("Basal", "Predators")){
      if (length(pred_labels) > 4) {
        x = c(1.6, 0.3, 2.9, 0.3, 2.9, 1.6)
        y = c(1.6, 1.1, 1.1, 0.4, 0.4, 0.1)
      } else {
        x = c(0.3, 2.9, 0.3, 2.9, 1.6)
        y = c(1.1, 1.1, 0.4, 0.4, 0.1)
      }
    } else {
      x = c(1.6, 0.3, 2.9, 0.3, 2.9, 1.6)
      y = c(1.6, 1.1, 1.1, 0.4, 0.4, 0.1)
    }
  } else {
    x = c(0.8, 2.4, 0.3, 2.9, 0.3, 2.9, 1.6)
    y = c(1.6, 1.6, 1.1, 1.1, 0.4, 0.4, 0.1)
  }
  
  ly = matrix(c(x, y), ncol=2)
  
  # Generate the SEM plot
  sem.plot <- qgraph(EL,
                     edge.color = sem.df$color,
                     edge.labels = sem.df$Std.Estimate,
                     edge.label.position = 0.25,
                     edge.label.cex = edge.cex,
                     edge.label.margin = 0.01,
                     edge.label.bg = sem.df$color_bg,
                     edge.label.color = sem.df$color,
                     lty = sem.df$style,
                     layout = ly,
                     shape = "ellipse",
                     label.cex = 1.5,
                     label.scale = F,
                     # label.font = "TT Arial",
                     borders = F,
                     vsize = size[1],
                     vsize2 = size[2],
                     edge.width = ifelse(is.na(sem.df$Std.Estimate), 1, abs(sem.df$Std.Estimate) * 10),
                     mar = c(2,4,2,4), 
                     DoNotPlot = T
  )
  
  # Add response variable labels
  sem.labels <-data.frame(Reponse = res_labels,
                          R2 = "")
  sem.labels <-  rbind(sem.labels, data.frame(Reponse = pred_labels,
                                              R2 = paste("R2: ", summary(sem)$R2$R.squared)
  ))
  
  # Append labels to plot
  sem.plot$graphAttributes$Nodes$labels <- paste(sem.labels[, 1], 
                                                 sem.labels[, 2], 
                                                 sep = "\n")
  return(sem.plot)
}

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
# Perform Principal Component Analysis (PCA) on quality-related covariates
PCA.urban = prcomp(PCA.urban.df, center = T, scale. = T)
# Store the first principal component as the 'Urban' metric
fw.df$Urbanization <- scale(predict(PCA.urban)[, 1])

# Compute the proportion of predator species relative to total species richness
fw.df$Predator.frac <- fw.df$Predator / fw.df$Sp.Richness

# Run and summarize SEM models
Urban.SEM.decoupl = psem(lm(Predator.frac ~ Urbanization, fw.df),
                         lm(Frac.Cluster.BG ~ Urbanization, fw.df),
                         lm(Link.Density ~ Predator.frac + Urbanization, fw.df), 
                         lm(Mean.ShortestPath ~ Urbanization + Frac.Cluster.BG + Link.Density, fw.df),
                         lm(Modularity ~ Urbanization + Mean.ShortestPath + Predator.frac + Frac.Cluster.BG + Link.Density, fw.df))

summary(Urban.SEM.decoupl, .progressBar = F)

# Plot SEM models
Urban.SEM.decoupl.plot <- plot_sem(Urban.SEM.decoupl, res_labels = "Urbanization", pred_labels = c("Proportion\nof predators",
                                                                                                   "Aquatic-terrestrial\nmodule size",
                                                                                                   "Link density",
                                                                                                   "Mean shortest\npath",
                                                                                                   "Modularity"), 
                                   size = c(25, 18), edge.cex = 1.8)


plot(Urban.SEM.decoupl.plot)
