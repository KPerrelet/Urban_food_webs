path <- "C:/Users/perrelki/switchdrive/Institution/"
# path <- "C:/Users/Kilia/switchdrive2/Institution/"

setwd(paste0(path, "/Chapter4/"))

################################################################################
######################### food web metrics calculation #########################
################################################################################

library(dplyr)
library(plyr)
library(bipartite)
library(moments)
library(igraph)

# Set a random seed for reproducibility
set.seed(1)

# Function to calculate unweighted trophic level (TL), omnivory index, and coherence index
# Based on prey-averaged approach from Williams & Martinez (2004)
# Original code from Blackman et al., 2022
Unw.TL.Omn.q = function(DietMatrix, Spe.Richness. = Spe.Richness, basal. = basal){
  TL1 = -t(DietMatrix)/apply(DietMatrix, 2, sum) # Normalize diet matrix
  TL1[basal.,] = 0 # Set basal species TL to 0
  diag(TL1) = 1 # Set diagonal elements to 1 (self-dependency)
  TL2 = c(rep(1, Spe.Richness.)) # Vector for solving linear equations
  if(det(TL1) < 1e-14) {TL = NA; Omn = NA; q = NA # Avoid error from computational singular matrix.
  } else {
    TL = solve(TL1, TL2)

    # Calculate omnivory index (variation in TL across prey)
    Omn = rep(NA, Spe.Richness.)
    for(s in 1:Spe.Richness.){
      if(s %in% basal.) next
      if(apply(DietMatrix, 2, sum)[s] == 1){ Omn[s] = 0
      }else{
        Omn[s] = sd(TL[which(DietMatrix[, s]==1)])
      }
    }

    # Calculate coherence index (trophic level consistency)
    Distance1 = matrix(TL, Spe.Richness., Spe.Richness., byrow = T)
    Distance2 = t(Distance1)
    Distance = DietMatrix * (Distance2 - Distance1)
    q = sd(Distance[which(Distance != 0)])
  } # end of else.

  return(list(TL = TL, Omn = Omn, q = q))
}

# Load data
taxamat.df <- read.csv("taxamat.csv")
taxo.df <- read.csv("taxo.csv")
MW <- read.csv("metaweb_processed.csv")

# Define a list of resource categories based on resources available in the metaweb
resources <- c("Aquatic detritus", "Terrestrial detritus", "Plants", "Algae", "Plankton", "Fungi", "Microbes", "Scavenger", "Haematophagous", "Coprophagous")
for (res in resources) {
  taxamat.df[, ncol(taxamat.df) + 1] <- 1
  names(taxamat.df)[ncol(taxamat.df)] <- res
  taxo.df <- rbind.fill(taxo.df, data.frame(taxa = res,
                                            Family = res))
  
}

# Aggregate the taxonomic matrix by site (right now, 1 site = aquatic + terrestrial communities)
bg.df <- aggregate(taxamat.df[, -which(colnames(taxamat.df) %in% c("site_no", "type"))], by = list(taxamat.df$site_no), FUN = sum)
colnames(bg.df)[1] <- "site_no"

# Create an adjacency matrix for the metaweb
MW.mat <- matrix(0,
                 nrow = length(unique(c(MW$Target_Name, MW$Source_Name))),
                 ncol = length(unique(c(MW$Target_Name, MW$Source_Name))),
                 dimnames = list(sort(unique(c(MW$Target_Name, MW$Source_Name))), sort(unique(c(MW$Target_Name, MW$Source_Name)))))
MW.mat[rownames(table(MW$Target_Name, MW$Source_Name)),
       colnames(table(MW$Target_Name, MW$Source_Name))] <- table(MW$Target_Name, MW$Source_Name)
MW.mat <- ifelse(MW.mat > 0, 1, 0) # Convert counts to binary presence/absence

# Iterate through each site to calculate food web metrics
res <- data.frame()
for(i in 1:nrow(bg.df)){ 
  site_no <- bg.df[i, ] # Extract site information
  sp <- colnames(bg.df)[-1][which(bg.df[i, 2:length(bg.df)] >= 1)] # Identify present species
  
  MW.local <- MW.mat[which(rownames(MW.mat) %in% sp),
                     which(colnames(MW.mat) %in% sp)] # Subset interaction matrix
  
  sp <- sp[which(sp %in% colnames(MW.mat))] # Ensure species exist in the metaweb
  
  # Identify isolated nodes (species with no interactions)
  Iso.Node <- names(which(apply(MW.local, 1, sum) == 0 & apply(MW.local, 2, sum) == 0))

  if (length(Iso.Node) != 0) {
    sp = setdiff(sp, Iso.Node) # Remove isolated nodes from species list
    MW.local <- MW.mat[which(rownames(MW.mat) %in% sp), 
                       which(colnames(MW.mat) %in% sp)] # Update interaction matrix
  }

  # Compute key food-web metrics
  Sp.Richness = length(sp[!(sp %in% resources)])  # Exclude resources

  Skew.Degree <- moments::skewness(igraph::degree(graph.adjacency(MW.local, mode = "directed")))
  No.Link = sum(MW.local)
  Link.Density = No.Link/Sp.Richness
  Connectance = No.Link/(Sp.Richness^2)
  Nic.over.Con.Jac = networklevel(MW.local, index = "niche overlap", dist = "jaccard")[1] # viewed as consumers vs. resources but actually the same set of taxa.
  Mean.Gen = mean(apply(MW.local, 2, sum))
  
  # Trophic level relevant metrics
  TL.info = Unw.TL.Omn.q(MW.local, ncol(MW.local), which(apply(MW.local, 2, sum)==0))
  Mean.TL = mean(TL.info$TL)
  Omnivory = mean(TL.info$Omn, na.rm = T) # In terms of TL, not primary producer
  Coherence = TL.info$q # Lower q, higher trophic coherence.
  Basal = length(TL.info$TL[TL.info$TL <= 2 & TL.info$TL > 1])
  Predator = length(TL.info$TL[TL.info$TL > 2])
  
  # Topological metrics
  Nestedness <- unname(nestednodf(MW.local)$statistic)[3]

  # Initialize storing vectors
  Size.Cluster.BG <- c()
  Frac.Cluster.BG <- c()
  Modularity <- c()
  for (j in c(1:100)) { # 100 iterations as the cluster_louvain algorithm is heuristic
    g <- graph_from_adjacency_matrix(MW.local, mode = "undirected") 
    
    # Determine community structure
    Community_rep <- cluster_louvain(g) 
    
    # Calculate modularity and append it to the storing vector
    Modularity <- append(Modularity, modularity(Community_rep)) 

    # Identify the clusters containing terrestrial and aquatic detritus, respectively
    ID.Cluster.BG <- which(membership(Community_rep) == membership(Community_rep)["Terrestrial detritus"] |
                             membership(Community_rep) == membership(Community_rep)["Aquatic detritus"]) 

    # Calculate the size and fraction of these clusters
    Size.Cluster.BG <- append(Size.Cluster.BG,
                              length(ID.Cluster.BG) - length(which(names(ID.Cluster.BG) %in% resources))) 
    Frac.Cluster.BG <- append(Frac.Cluster.BG,
                              (length(ID.Cluster.BG) - length(which(names(ID.Cluster.BG) %in% resources))) / Sp.Richness)
    }

  Modularity <- mean(Modularity)
  Size.Cluster.BG <- mean(Size.Cluster.BG)
  Frac.Cluster.BG <- mean(Frac.Cluster.BG)
  
  # Calculate shortest path lengths
  # We remove resources to avoid very distance species being connected by going through broadly defined and shared resources (e.g., plants)
  paths <- distances(graph_from_adjacency_matrix(MW.local[-which(rownames(MW.local) %in% resources),
                                                          -which(colnames(MW.local) %in% resources)], 
                                                 mode = "undirected"),
                     mode = "out")
  # For the species that do not connect, we take the maximum and add 1 (to simulate a connection through a shared resource)
  paths[is.infinite(paths) == T] <- max(paths[is.infinite(paths) == F]) + 1
  Mean.ShortestPath <- mean(paths)

  ### Null model
  
  # Create empty matrix the size of the local matrix
  MW.local.rand <-  matrix(0,
                             ncol = ncol(MW.local),
                             nrow = nrow(MW.local))

  # Randomly add 1 in the matrix, except for resources (which cannot eat other species)
  MW.local.rand[sample(which(!(col(MW.local.rand) %in% which(colnames(MW.local) %in% resources))),
                         sum(MW.local),
                         replace = F)] <- 1
  colnames(MW.local.rand) <- colnames(MW.local)
  rownames(MW.local.rand) <- rownames(MW.local)

  # Calculate the same metrics as before
  Skew.Degree.rand <- skewness(igraph::degree(graph.adjacency(MW.local.rand, mode = "directed")))
  Nic.over.Con.Jac.rand = networklevel(MW.local.rand, index = "niche overlap", dist = "jaccard")[1] # viewed as consumers vs. resources but actually the same set of taxa.
  Mean.Gen.rand = mean(apply(MW.local.rand, 2, sum))
  
  TL.info.rand = Unw.TL.Omn.q(MW.local.rand, ncol(MW.local.rand), which(apply(MW.local.rand, 2, sum)==0))
  Mean.TL.rand = mean(TL.info.rand$TL)
  Omnivory.rand = mean(TL.info.rand$Omn, na.rm = T) # In terms of TL, not primary producer
  Coherence.rand = TL.info.rand$q # Lower q, higher trophic coherence.
  Basal.rand = length(TL.info.rand$TL[TL.info.rand$TL <= 2 & TL.info.rand$TL > 1])
  Predator.rand = length(TL.info.rand$TL[TL.info.rand$TL > 2])
  
  Nestedness.rand = unname(nestednodf(MW.local.rand)$statistic)[3]
  
  Size.Cluster.BG.rand <- c()
  Frac.Cluster.BG.rand <- c()
  Modularity.rand <- c()
  for (i in c(1:100)) {
    Community.rand <- cluster_louvain(graph_from_adjacency_matrix(MW.local.rand, mode = "undirected"))
    Modularity.rand <- append(Modularity.rand, modularity(Community.rand))

    ID.Cluster.BG.rand <- which(membership(Community.rand) == membership(Community.rand)["Terrestrial detritus"] |
                             membership(Community.rand) == membership(Community.rand)["Aquatic detritus"] 
    )
    
    Size.Cluster.BG.rand <- append(Size.Cluster.BG.rand, 
                                     length(ID.Cluster.BG.rand) - length(which(names(ID.Cluster.BG.rand) %in% resources)))
    Frac.Cluster.BG.rand <- append(Frac.Cluster.BG.rand, 
                                     (length(ID.Cluster.BG.rand) - length(which(names(ID.Cluster.BG.rand) %in% resources))) / Sp.Richness)
  }
  
  Modularity.rand <- mean(Modularity.rand)
  Size.Cluster.BG.rand <- mean(Size.Cluster.BG.rand)
  Frac.Cluster.BG.rand <- mean(Frac.Cluster.BG.rand)
  
  paths <- distances(graph_from_adjacency_matrix(MW.local.rand[-which(rownames(MW.local.rand) %in% resources),
                                                          -which(colnames(MW.local.rand) %in% resources)], 
                                                 mode = "undirected"),
                     mode = "out")
  paths[is.infinite(paths) == T] <- max(paths[is.infinite(paths) == F]) + 1
  Mean.ShortestPath.rand <- mean(paths)
  
  # Fill in dataframe
  temp.data = data.frame(
    site_no = site_no$site_no,
    No.Iso.Node = length(Iso.Node),
    # No.NoRes.Con = length(NoRes.Con),
    Sp.Richness = Sp.Richness,
    Basal = Basal,
    Predator = Predator, 
    No.Link = No.Link,
    Link.Density = Link.Density,
    Connectance = Connectance,
    Nic.over.Con.Jac = Nic.over.Con.Jac,
    Mean.Gen = Mean.Gen,
    Mean.TL = Mean.TL,
    Omnivory = Omnivory,
    Coherence = Coherence,
    Skew.Degree = Skew.Degree, 
    Nestedness = Nestedness,
    Modularity = Modularity,
    Size.Cluster.BG = Size.Cluster.BG,
    Frac.Cluster.BG = Frac.Cluster.BG,
    Mean.ShortestPath = Mean.ShortestPath,
    
    # Null model metrics 
    Basal.rand = Basal.rand,
    Predator.rand = Predator.rand,
    Skew.Degree.rand = Skew.Degree.rand,
    Nic.over.Con.Jac.rand = Nic.over.Con.Jac.rand,
    Mean.Gen.rand = Mean.Gen.rand,
    Mean.TL.rand = Mean.TL.rand,
    Omnivory.rand = Omnivory.rand,
    Coherence.rand = Coherence.rand,
    Nestedness.rand = Nestedness.rand,
    Modularity.rand = Modularity.rand,
    Size.Cluster.BG.rand = Size.Cluster.BG.rand, 
    Frac.Cluster.BG.rand = Frac.Cluster.BG.rand,
    Mean.ShortestPath.rand = Mean.ShortestPath.rand
  )
  res <- rbind(res, temp.data)
}

write.csv(res, "foodweb_metrics.csv", row.names = F)
