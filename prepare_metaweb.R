# path <- "C:/Users/perrelki/switchdrive/Institution/"
path <- "C:/Users/Kilia/switchdrive2/Institution/"

setwd(paste0(path, "/Chapter4/"))

minsize <- 10
scheme <- "strict_"
# scheme <- ""


################################################################################
############################## Metaweb preparation #############################
################################################################################
library(dplyr)
library(plyr)

# Load data
taxo.df <- read.csv("taxo.csv")
MW <- read.csv("metaweb.csv")

# Add resources to taxonomic dataset
resources <- c("Aquatic detritus", "Terrestrial detritus", "Plants", "Algae", "Plankton", "Fungi", "Microbes", "Scavenger", "Haematophagous", "Coprophagous")
for (res in resources) {
  taxo.df <- rbind.fill(taxo.df, data.frame(taxa = res,
                                            Family = res))
  
}

# Add a specific interaction entry manually to the metaweb dataset
MW <- rbind(MW, data.frame(Interaction = "Caeciliusidae Terrestrial detritus",
                           Source_Name = "Caeciliusidae",
                           Target_Name = "Terrestrial detritus",
                           Source_Rank = "Family",
                           Target_Rank = "Feeding guild",
                           ID = 100,
                           Source_Life_Stage = "Young and Adult",
                           Target_Life_Stage = NA))

# Expand the metaweb dataset by resolving taxonomic ambiguities for specific taxa
for (taxa in c("Collembola", "Annelida", "Acari", "Isopoda")) {
  for (source in MW[MW$Target_Name == taxa, ]$Source_Name) {
    if (taxa == "Collembola") {
      for (target in na.omit(unique(taxo.df[taxo.df$Class == taxa, ]$Family))) {
        MW <- rbind(MW, data.frame(Interaction = paste(source, target, sep = " "),
                                   Source_Name = source,
                                   Target_Name = target,
                                   Source_Rank = "Family",
                                   Target_Rank = "Family",
                                   ID = NA,
                                   Source_Life_Stage = MW[MW$Target_Name == taxa & MW$Source_Name == source, ]$Source_Life_Stage,
                                   Target_Life_Stage = NA))
      }
      MW <- MW[-which(MW$Target_Name == taxa & MW$Source_Name == source), ] 
    }
    if (taxa == "Annelida") {
      for (target in na.omit(unique(taxo.df[taxo.df$Phylum == taxa, ]$Family))) {
        MW <- rbind(MW, data.frame(Interaction = paste(source, target, sep = " "),
                                   Source_Name = source,
                                   Target_Name = target,
                                   Source_Rank = "Family",
                                   Target_Rank = "Family",
                                   ID = NA,
                                   Source_Life_Stage = MW[MW$Target_Name == taxa & MW$Source_Name == source, ]$Source_Life_Stage,
                                   Target_Life_Stage = NA))
      }
      # Remove original unresolved interactions
      MW <- MW[-which(MW$Target_Name == taxa & MW$Source_Name == source), ] 
    }
    if (taxa == "Acari") {
      for (target in na.omit(unique(taxo.df[taxo.df$Order == "Mesostigmata", ]$Family))) {
        MW <- rbind(MW, data.frame(Interaction = paste(source, target, sep = " "),
                                   Source_Name = source,
                                   Target_Name = target,
                                   Source_Rank = "Family",
                                   Target_Rank = "Family",
                                   ID = NA,
                                   Source_Life_Stage = MW[MW$Target_Name == taxa & MW$Source_Name == source, ]$Source_Life_Stage,
                                   Target_Life_Stage = NA))
      }
      MW <- MW[-which(MW$Target_Name == taxa & MW$Source_Name == source), ] 
    }
    if (taxa == "Isopoda") {
      for (target in na.omit(unique(taxo.df[taxo.df$Order == "Isopoda", ]$Family))) {
        MW <- rbind(MW, data.frame(Interaction = paste(source, target, sep = " "),
                                   Source_Name = source,
                                   Target_Name = target,
                                   Source_Rank = "Family",
                                   Target_Rank = "Family",
                                   ID = NA,
                                   Source_Life_Stage = MW[MW$Target_Name == taxa & MW$Source_Name == source, ]$Source_Life_Stage,
                                   Target_Life_Stage = NA))
      }
      MW <- MW[-which(MW$Target_Name == taxa & MW$Source_Name == source), ] 
    }
  }
}

# Merge metaweb with taxonomy dataset based on different taxonomic ranks
MW <- MW %>%
  left_join(taxo.df, by = c("Source_Name" = "Phylum")) %>%
  left_join(taxo.df, by = c("Source_Name" = "Class")) %>%
  left_join(taxo.df, by = c("Source_Name" = "Order")) %>%
  dplyr::select(Interaction, Source_Name, Target_Name, Source_Rank, Target_Rank, Source_Life_Stage, Target_Life_Stage, Family.x, Family.y, Family)

# Replace higher taxonomic ranks with family names
MW[MW$Source_Name == "Annelida", ]$Source_Name <- MW[MW$Source_Name == "Annelida", ]$Family.x
MW[MW$Source_Name == "Collembola", ]$Source_Name <- MW[MW$Source_Name == "Collembola", ]$Family.y
MW[MW$Source_Name == "Isopoda", ]$Source_Name <- MW[MW$Source_Name == "Isopoda", ]$Family

# Remove duplicate rows
MW <- unique(MW)

# Assign functional groups (feeding guilds) based on target interactions
MW$fg <- NA
for (family in unique(MW$Source_Name)) {
  
  fg <- ifelse(MW[MW$Source_Name == family, ]$Target_Name %in% c(resources), 
               ifelse(MW[MW$Source_Name == family, ]$Target_Name == "Plants", "Herbivore", "Detritivore"), 
               "Predator")
  
  if (length(unique(fg)) > 1 & "Predator" %in% fg) {
    MW[MW$Source_Name == family, ]$fg <- "Omnivore"
  } else{
    if ("Detritivore" %in% fg & "Herbivore" %in% fg) {
      MW[MW$Source_Name == family, ]$fg <- "Phytosaprophage"
    } else {
      MW[MW$Source_Name == family, ]$fg <- unique(fg)
    }
  }
}

# Remove unnecessary columns
MW <- MW %>%
  select(-c("Family", "Family.x", "Family.y", "Source_Life_Stage", "Target_Life_Stage"))

# Save processed metaweb to a CSV file
write.csv(MW, "metaweb_processed.csv", row.names = F)
