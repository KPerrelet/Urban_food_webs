################################################################################
#################################### fig S4 ####################################
################################################################################

library(piecewiseSEM)
library(ggplot2)
library(ggpubr)
library(semEff)

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

# Get direct and indirect effects of the urbanization covariate (bootstrapping-based method, 10000 iterations)
Urban.eff.decoupl <- get_sem_eff(Urban.SEM.decoupl, "Urbanization", R = 10000)

Urban.eff.decoupl$effect_type <- factor(Urban.eff.decoupl$effect_type, levels=c("direct", "total", "indirect"))

# Create coefficient plot
ggplot(Urban.eff.decoupl[Urban.eff.decoupl$effect_type != "total", ], aes(y = response, x = effect, fill = predictor, color = predictor, shape = effect_type)) + 
  geom_tile(data = Urban.eff.decoupl[Urban.eff.decoupl$effect_type != "total" & Urban.eff.decoupl$response %in% c(
    "Predator.frac", 
    "Link.Density", 
    "Modularity"
  ), ], 
  aes(y = response, x = 0, height = 1, width = Inf), fill =  "grey95", color = NA) +  geom_vline(xintercept = 0) +
  geom_point(size = 3, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci), 
                width = 0, linewidth = 1.2, position = position_dodge(width=0.9)) + 
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(16, 17), labels = c("Direct", "Indirect")) + 
  scale_y_discrete(
    limits = c("Modularity", 
               "Mean.ShortestPath",
               "Link.Density",
               "Frac.Cluster.BG",
               "Predator.frac"),
    labels = c("Modularity", 
               "Mean shortest\ndistance", 
               "Link density",
               "Aquatic-terrestrial\nmodule relative sizes",
               "Proportion of predators")) +
  labs(x = "Standardized coefficients", y = "") + 
  guides(y = ggh4x::guide_axis_nested(delim = "&"), 
         fill = "none",
         color = "none", 
         shape = guide_legend(title = "Effect type")) + 
  theme_pubr(base_size = 18) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 18), 
        axis.title = element_text(size = 18), 
        ggh4x.axis.nesttext.y = element_text(angle = 90, hjust = .5, vjust = 1))  
