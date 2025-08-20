################################################################################
##################################### fig 4 #################################### 
################################################################################

library(piecewiseSEM)
library(semEff)
library(ggplot2)
library(ggpubr)
library(scales)

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

# Get direct and indirect effects of each land-use covariate (bootstrapping-based method, 10000 iterations)
Urban.eff.struct <- get_sem_eff(Urban.SEM.struct, "Urbanization", R = 10000)
Urban.eff.compo <- get_sem_eff(Urban.SEM.compo, "Urbanization", R = 10000)
QualityQuantity.eff.struct <- get_sem_eff(QuantityQuality.SEM.struct, c("Quantity", "Quality"), R = 10000)
QualityQuantity.eff.compo <- get_sem_eff(QuantityQuality.SEM.compo, c("Quantity", "Quality"), R = 10000)
ConnCompl.eff.struct <- get_sem_eff(ConnCompl.SEM.struct, c("Connectivity", "Complementarity"), R = 10000)
ConnCompl.eff.compo <- get_sem_eff(ConnCompl.SEM.compo, c("Connectivity", "Complementarity"), R = 10000)

# Store all direct and indirect effects into a dataframe
sem.eff <- do.call("rbind", list(Urban.eff.struct, 
                                 Urban.eff.compo, 
                                 QualityQuantity.eff.struct, 
                                 QualityQuantity.eff.compo, 
                                 ConnCompl.eff.struct, 
                                 ConnCompl.eff.compo))

# Reorder the dataframe
sem.eff$predictor <- factor(sem.eff$predictor, levels=c("Urbanization", "Quantity", "Quality", "Connectivity", "Complementarity"))

# Create two groups for visualization purposes
sem.eff$group <- ifelse(sem.eff$response %in% c("Sp.Richness", "Connectance", "Modularity", "Nestedness", "Nic.over.Con.Jac"),
                        "Structure", "Composition")

sem.eff$response <- paste(sem.eff$response, sem.eff$group, sep = "&")

sem.eff$effect_type <- factor(sem.eff$effect_type, levels=c("direct", "indirect"))

# Create coefficient plot
ggplot(sem.eff, aes(y = response, x = effect, fill = predictor, color = predictor, shape = effect_type)) + 
  geom_tile(data = sem.eff[sem.eff$response %in% c(
    "Nestedness&Structure",
    "Connectance&Structure",
    "Coherence&Composition", 
    "Mean.Gen&Composition",
    "Skew.Degree&Composition"), ], 
    aes(y = response, x = 0, height = 1, width = Inf), fill =  "grey95", color = NA) +  geom_vline(xintercept = 0) +
  geom_point(aes(alpha = ifelse(lower_ci <= 0 & upper_ci >= 0, "insign", "sign")), size = 4, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci, alpha = ifelse(lower_ci <= 0 & upper_ci >= 0, "insign", "sign")), 
                width = 0, linewidth = 1.2, position = position_dodge(width=0.9)) + 
  facet_wrap(~ predictor, nrow = 1, scales = "free_x") + 
  scale_color_manual(values = c("#636363", "#a6d96a", "#1b7837", "#92c5de", "#2166ac")) +
  scale_fill_manual(values = c("#636363", "#a6d96a", "#1b7837", "#92c5de", "#2166ac")) +
  scale_alpha_manual(values = c(0.3, 1)) +
  # scale_shape_manual(values = c(16, 15, 17), labels = c("Direct", "Total", "Indirect")) + 
  scale_shape_manual(values = c(16, 17), labels = c("Direct", "Indirect")) + 
  scale_x_continuous(labels = label_number(accuracy = 0.1)) + 
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
               "Trophic\nlevel&Composition",
               "Degree\nskewness&Composition")) +
  labs(x = "Standardized coefficients and 95% CI", y = "") + 
  guides(y = ggh4x::guide_axis_nested(delim = "&"), 
         fill = "none",
         color = "none", 
         alpha = "none",
         shape = guide_legend(title = "Effect type")) + 
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), 
        strip.text = element_text(size = 18), 
        ggh4x.axis.nesttext.y = element_text(angle = 90, hjust = .5, vjust = 1))  
