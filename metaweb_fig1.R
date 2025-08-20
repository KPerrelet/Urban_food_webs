################################################################################
################################# Graph summary ################################
################################################################################
library(dplyr)
library(plyr)
library(igraph)
library(ggplot2)

# Load data
MW <- read.csv("metaweb_processed.csv")
taxo.df <- read.csv("taxo.csv")

# Define a list of resource categories to be included in the dataset
resources <- c("Aquatic detritus", "Terrestrial detritus", "Plants", "Algae", "Plankton", "Fungi", "Microbes", "Scavenger", "Haematophagous", "Coprophagous")
# Append these resources to the taxonomic data frame
for (res in resources) {
  taxo.df <- rbind.fill(taxo.df, data.frame(taxa = res,
                                            Family = res))
}

# Merge taxonomic information with the metaweb dataset based on source names
MW <- merge(MW, unique(taxo.df[, c(2:5)]), by.x = "Source_Name", by.y = "Family")
colnames(MW)[8:10] <- paste(colnames(MW)[8:10], "Source", sep = "_")
# Merge taxonomic information with the metaweb dataset based on target names
MW <- merge(MW, unique(taxo.df[, c(2:5)]), by.x = "Target_Name", by.y = "Family")
colnames(MW)[11:13] <- paste(colnames(MW)[11:13], "Target", sep = "_")
# Fill missing taxonomic classifications with available target names
MW[is.na(MW$Phylum_Target), ]$Phylum_Target <- MW[is.na(MW$Phylum_Target), ]$Target_Name
MW[is.na(MW$Order_Target), ]$Order_Target <- MW[is.na(MW$Order_Target), ]$Target_Name
MW[is.na(MW$Class_Target), ]$Class_Target <- MW[is.na(MW$Class_Target), ]$Target_Name

# Assign taxonomic group names based on different classification levels
MW$Grp_Name_Source <- MW$Order_Source
MW[-which(MW$Class_Source == "Insecta"), ]$Grp_Name_Source <- MW[-which(MW$Class_Source == "Insecta"), ]$Class_Source
MW[which(MW$Phylum_Source != "Arthropoda"), ]$Grp_Name_Source <- MW[which(MW$Phylum_Source != "Arthropoda"), ]$Phylum_Source

# Special classification rules for Arachnida, Crustacea, and Myriapoda
MW[MW$Class_Source == "Arachnida" | MW$Class_Source == "Opiliones", ]$Grp_Name_Source <- ifelse(MW[MW$Class_Source == "Arachnida" | MW$Class_Source == "Opiliones", ]$Order_Source == "Araneae", "Araneae", "Acari")
MW[MW$Grp_Name_Source == "Hexanauplia" | 
     MW$Grp_Name_Source == "Branchiopoda" | 
     MW$Grp_Name_Source == "Malacostraca" |
     MW$Grp_Name_Source == "Ostracoda", ]$Grp_Name_Source <- "Crustacea"
MW[MW$Grp_Name_Source == "Chilopoda" | 
     MW$Grp_Name_Source == "Symphyla" | 
     MW$Grp_Name_Source == "Diplopoda", ]$Grp_Name_Source <- "Myriapoda"
MW[MW$Grp_Name_Source == "Megaloptera" | 
     MW$Grp_Name_Source == "Neuroptera", ]$Grp_Name_Source <- "Neuropterida"

# Repeat similar processing for target groups
MW$Grp_Name_Target <- MW$Order_Target
MW[-which(MW$Class_Target == "Insecta"), ]$Grp_Name_Target <- MW[-which(MW$Class_Target == "Insecta"), ]$Class_Target
MW[which(MW$Phylum_Target != "Arthropoda"), ]$Grp_Name_Target <- MW[which(MW$Phylum_Target != "Arthropoda"), ]$Phylum_Target
MW[MW$Class_Target == "Arachnida" | MW$Class_Target == "Opiliones", ]$Grp_Name_Target <- ifelse(MW[MW$Class_Target == "Arachnida" | MW$Class_Target == "Opiliones", ]$Order_Target == "Araneae", "Araneae", "Acari")
MW[MW$Grp_Name_Target == "Hexanauplia" | 
     MW$Grp_Name_Target == "Branchiopoda" | 
     MW$Grp_Name_Target == "Malacostraca" |
     MW$Grp_Name_Target == "Ostracoda", ]$Grp_Name_Target <- "Crustacea"
MW[MW$Grp_Name_Target == "Chilopoda" | 
     MW$Grp_Name_Target == "Symphyla" | 
     MW$Grp_Name_Target == "Diplopoda", ]$Grp_Name_Target <- "Myriapoda"
MW[MW$Grp_Name_Target == "Megaloptera" | 
     MW$Grp_Name_Target == "Neuroptera", ]$Grp_Name_Target <- "Neuropterida"

MW[which(MW$Grp_Name_Target %in% c("Algae",
                                   "Aquatic detritus", 
                                   "Coprophagous",
                                   "Haematophagous",
                                   "Fungi", 
                                   "Microbes", 
                                   "Plankton", 
                                   "Plants",
                                   "Scavenger", 
                                   "Terrestrial detritus")), ]$Grp_Name_Target <- "Basal"

# Extract relevant columns for constructing adjacency matrix
MW.sub <- MW[, which(colnames(MW) %in% c("Grp_Name_Source", "Grp_Name_Target", "fg"))]

# Create adjacency matrix representing interactions between groups
MW.full.mat <- matrix(0, 
                      nrow = length(unique(c(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source))), 
                      ncol = length(unique(c(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source))),
                      dimnames = list(sort(unique(c(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source))), 
                                      sort(unique(c(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source)))))

MW.full.mat[rownames(table(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source)), 
            colnames(table(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source))] <- table(MW.sub$Grp_Name_Target, MW.sub$Grp_Name_Source)

# Convert adjacency matrix into a graph object
full.fw <- graph.adjacency(MW.full.mat, mode = "directed")

# Create layout and node properties for visualization
fw.layout.nodes <- as.data.frame(layout_in_circle(full.fw))
fw.layout.nodes$Grp <- colnames(MW.full.mat)
fw.layout.nodes$fg_percent <- NA
fw.layout.nodes$Phylum <- NA
fw.layout.nodes$size <- NA
for (i in 1:nrow(fw.layout.nodes)) {
  node <- fw.layout.nodes[i, ]
  if (node$Grp != "Basal") {
    fw.layout.nodes[i, ]$size <- length(unique(MW[MW$Grp_Name_Source == node$Grp, ]$Source_Name))
    fw.layout.nodes[i, ]$Phylum <- unique(MW[MW$Grp_Name_Source == node$Grp, ]$Phylum_Source)
    fw.layout.nodes[i, ]$fg_percent <- nrow(MW[MW$Grp_Name_Source == node$Grp & MW$fg %in% c("Omnivore"), ]) / nrow(MW[MW$Grp_Name_Source == node$Grp, ])
  } else {
    fw.layout.nodes[i, ]$Phylum <- "Basal"
    fw.layout.nodes[i, ]$fg_percent <- 0
    fw.layout.nodes[i, ]$size <- 1
  }
}

# Define edge properties for visualization
fw.layout.edges <- unique(as.data.frame(as_edgelist(full.fw)))
colnames(fw.layout.edges) <- c("from", "to")
fw.layout.edges <- merge(fw.layout.edges, fw.layout.nodes, by.x = "from", by.y = "Grp")
fw.layout.edges <- merge(fw.layout.edges, fw.layout.nodes, by.x = "to", by.y = "Grp", suffixes = c(".from", ".to"))
fw.layout.edges$thickness <- NA
for (i in 1:nrow(fw.layout.edges)) {
  edge <- fw.layout.edges[i, ]
  fw.layout.edges[i, ]$thickness <- MW.full.mat[edge$from, 
                                                edge$to] 
}

fw.layout.edges <- fw.layout.edges[-which(fw.layout.edges$from == fw.layout.edges$to), ]

fw.layout.nodes$position_x = c(rep("left", 6), rep("right", 12), rep("left", 5))
fw.layout.nodes$position_y = c(rep("top", 12), rep("bottom", 11))

# Create food web visualization 
ggplot() +
  geom_segment(data = fw.layout.edges, aes(x = V1.from, y = V2.from, xend = V1.to, yend = V2.to, linewidth = thickness), color = "grey") +
  geom_point(data=fw.layout.nodes, aes(x = V1, y = V2, color = Phylum, fill = fg_percent, size = size),  shape=21, stroke = 3) +
  geom_label(data=fw.layout.nodes[fw.layout.nodes$position_x == "left" & fw.layout.nodes$position_y == "top", ],
                   aes(x = V1*1.1, y = V2*1.1, label = Grp), size = 12, nudge_x = 0.2) +
  geom_label(data=fw.layout.nodes[fw.layout.nodes$position_x == "right" & fw.layout.nodes$position_y == "top", ],
                   aes(x = V1*1.1, y = V2*1.1, label = Grp), size = 12, nudge_x = -0.2) +
  geom_label(data=fw.layout.nodes[fw.layout.nodes$position_x == "left" & fw.layout.nodes$position_y == "bottom", ],
                   aes(x = V1*1.1, y = V2*1.1, label = Grp), size = 12, nudge_x = 0.2) +
  geom_label(data=fw.layout.nodes[fw.layout.nodes$position_x == "right" & fw.layout.nodes$position_y == "bottom", ],
            aes(x = V1*1.1, y = V2*1.1, label = Grp), size = 12, nudge_x = -0.2) +
    scale_x_continuous(expand=c(0, 0.3)) +
  scale_y_continuous(expand=c(0, 0.3)) + 
  scale_linewidth(range = c(1,7)) + 
  scale_size(range = c(5,30), 
             breaks = c(1, 40)) +
  scale_fill_distiller(type = "seq",
                        direction = 1,
                        palette = "Greys") +
  theme_void() + 
  theme(text = element_text(size = 32), 
        legend.key.size = unit(1, 'cm'), 
        legend.margin = margin(50),
        legend.spacing.y = unit(0.5, 'cm'), 
        legend.spacing.x = unit(0.5, 'cm')) + 
  guides(color = guide_legend(title = "Phylum", byrow = T, title.vjust = 2, override.aes = list(size = 8), order = 4), 
         fill = guide_legend(title = "Omnivory levels", byrow = T, title.vjust = 2, override.aes = list(size = 8), order = 3), 
         linewidth = guide_legend(title = "Number of links", byrow = F, title.vjust = 1, order = 2), 
         size = guide_legend(title = "Number of families", order = 1))


