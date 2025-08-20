################################################################################
##################################### fig 4 #################################### 
################################################################################

library(piecewiseSEM)
library(igraph)
library(qgraph)
library(scales)

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

# Function to extract SEM effects
get_sem_eff <- function(sem, predictors, R = 10000, seed = 123, parallel = "multicore", type = "nonparametric") {
  sem.boot <- bootEff(sem, R = R, seed = seed, parallel = parallel, type = type)
  sem.eff <- semEff(sem.boot, predictor = predictors, ci.type = "norm")$Effects$Table
  
  # Filter results to include only direct and indirect effects
  sem.eff <- 
    sem.eff[sem.eff$predictor %in% predictors & sem.eff$effect_type %in% c("direct", "indirect"), ]
  return(sem.eff)
}

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
fw.df$Urbanization <- scale(predict(PCA.urban)[, 1])

# Run and summarize SEM models
Urban.SEM.compo <- psem(lm(Skew.Degree ~ Urbanization, fw.df), 
                        lm(Mean.TL ~ Skew.Degree + Urbanization, fw.df), 
                        lm(Mean.Gen ~ Skew.Degree + Mean.TL + Urbanization, fw.df),
                        lm(Omnivory ~ Skew.Degree + Mean.TL + Urbanization, fw.df),
                        lm(Coherence ~ Skew.Degree + Mean.TL + Omnivory + Urbanization, fw.df))
summary(Urban.SEM.compo)

QuantityQuality.SEM.compo <-  psem(lm(Skew.Degree ~ Quantity + Quality, fw.df),
                                   lm(Mean.TL ~ Skew.Degree + Quantity + Quality, fw.df),
                                   lm(Mean.Gen ~ Skew.Degree + Mean.TL + Quantity + Quality, fw.df),
                                   lm(Omnivory ~ Skew.Degree + Mean.TL + Quantity + Quality, fw.df),
                                   lm(Coherence ~ Skew.Degree + Mean.TL + Omnivory + Quantity + Quality, fw.df))
summary(QuantityQuality.SEM.compo)

ConnCompl.SEM.compo =  psem(lm(Skew.Degree ~ Connectivity + Complementarity, fw.df),
                            lm(Mean.TL ~ Skew.Degree + Connectivity + Complementarity, fw.df),
                            lm(Mean.Gen ~ Skew.Degree + Mean.TL + Connectivity + Complementarity, fw.df),
                            lm(Omnivory ~ Skew.Degree + Mean.TL + Connectivity + Complementarity, fw.df),
                            lm(Coherence ~ Skew.Degree + Mean.TL + Omnivory + Connectivity + Complementarity, fw.df))
summary(ConnCompl.SEM.compo)

Urban.SEM.struct = psem(glm(Sp.Richness ~ Urbanization, family = "poisson", fw.df),
                        lm(Connectance ~ Urbanization + Sp.Richness, fw.df),
                        lm(Modularity ~ Urbanization + Sp.Richness + Connectance, fw.df),
                        lm(Nestedness ~ Urbanization + Sp.Richness + Connectance, fw.df),
                        lm(Nic.over.Con.Jac ~ Nestedness + Modularity + Urbanization + Sp.Richness + Connectance, fw.df)
)
summary(Urban.SEM.struct)
QuantityQuality.SEM.struct = psem(glm(Sp.Richness ~ Quantity + Quality, family = "poisson", fw.df),
                                  lm(Connectance ~ Quantity + Quality + Sp.Richness, fw.df),  
                                  lm(Modularity ~ Quantity + Quality + Sp.Richness + Connectance, fw.df),
                                  lm(Nestedness ~ Quantity + Quality + Sp.Richness + Connectance, fw.df),
                                  lm(Nic.over.Con.Jac ~ Nestedness + Modularity + Quantity + Quality + Sp.Richness + Connectance, fw.df)
)
summary(QuantityQuality.SEM.struct)

ConnCompl.SEM.struct = psem(glm(Sp.Richness ~ Connectivity + Complementarity, family = "poisson", fw.df),
                            lm(Connectance ~ Connectivity + Complementarity + Sp.Richness, fw.df),
                            lm(Modularity ~ Connectivity + Complementarity + Sp.Richness + Connectance, fw.df),
                            lm(Nestedness ~ Connectivity + Complementarity + Sp.Richness + Connectance, fw.df),
                            lm(Nic.over.Con.Jac ~ Nestedness + Modularity + Connectivity + Complementarity + Sp.Richness + Connectance, fw.df)
)
summary(ConnCompl.SEM.struct)

Urban.SEM.compo.plot <- plot_sem(Urban.SEM.compo, res_labels = "Urbanization", pred_labels = c("Degree\nskewness",
                                                                                               "Trophic\nlevel mean",
                                                                                               "Omnivory",
                                                                                               "Mean\ngenerality",
                                                                                               "Incoherence"))
Urban.SEM.struct.plot <- plot_sem(Urban.SEM.struct, res_labels = "Urbanization", pred_labels = c("Family\nrichness",
                                                                                                 "Connectance",
                                                                                                 "Modularity",
                                                                                                 "Nestedness",
                                                                                                 "Niche\noverlap"))

QuantityQuality.SEM.compo.plot <- plot_sem(QuantityQuality.SEM.compo, res_labels = c("Quantity", "Quality"), pred_labels = c("Degree\nskewness",
                                                                                                                             "Trophic\nlevel mean",
                                                                                                                             "Omnivory",
                                                                                                                             "Mean\ngenerality",
                                                                                                                             "Incoherence"))
QuantityQuality.SEM.struct.plot <- plot_sem(QuantityQuality.SEM.struct, res_labels = c("Quantity", "Quality"), pred_labels = c("Family\nrichness",
                                                                                                                               "Connectance",
                                                                                                                               "Modularity",
                                                                                                                               "Nestedness",
                                                                                                                               "Niche\noverlap"))

ConnCompl.SEM.compo.plot <- plot_sem(ConnCompl.SEM.compo, res_labels = c("Connectivity", "Complementarity"), pred_labels = c("Degree\nskewness",
                                                                                                                             "Trophic\nlevel mean",
                                                                                                                             "Omnivory",
                                                                                                                             "Mean\ngenerality",
                                                                                                                             "Incoherence"))
ConnCompl.SEM.struct.plot <- plot_sem(ConnCompl.SEM.struct, res_labels = c("Connectivity", "Complementarity"), pred_labels = c("Family\nrichness",
                                                                                                                               "Connectance",
                                                                                                                               "Modularity",
                                                                                                                               "Nestedness",
                                                                                                                               "Niche\noverlap"))

# Plot SEM models
plot(Urban.SEM.compo.plot)
plot(QuantityQuality.SEM.compo.plot)
plot(ConnCompl.SEM.compo.plot)
plot(Urban.SEM.struct.plot)
plot(QuantityQuality.SEM.struct.plot)
plot(ConnCompl.SEM.struct.plot)
