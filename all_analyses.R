## load packages
library(Morpho)
library(rgl)
library(tidyverse)
library(svglite)
library(MCMCglmm)
library(ape)
library(dispRity)
library(geomorph)
library(ggtree)
library(agricolae)
library(patchwork)
# set the working directory to the folder where landmark files are located
setwd('/Users/harveyjeffrey/project_data/landmarks')

## Detecting all files that are *.fcsv
landmark_files <- list.files(pattern = "*.fcsv")


## Creating an empty list to store the landmarks for each scan
scans_list <- list()

## Loading all the landmarks using a for loop (here it replaces
## "one_landmark_file" by each element in "landmark_files" sequentially)
for(one_landmark_file in landmark_files) {
  ## Reading the landmarks for the one file
  my_landmarks <- read.fcsv(one_landmark_file)
  ## Counting how many elements are present in the list
  number_of_elements <- length(scans_list)
  ## Adding the landmarks to the end of the list (+1)
  scans_list[[number_of_elements + 1]] <- my_landmarks
}

## Getting the name of each landmark file (i.e. removing the *.fcsv part)
landmark_files_names <- unlist(strsplit(landmark_files, split = "*.fcsv"))

## Adding the scan names to the list of landmarks
names(scans_list) <- landmark_files_names

## Converting the list into an array (an array is just a different way to store
## the data in R and is required by Morpho)
scans_array <- list2array(scans_list)

# state the landmark pairs
landmark_pairs <- as.matrix(read.table("landmark_pairs.txt"))

## Running the GPA and the PCA
procrustes_data <- procSym(scans_array, 
                          ## The pairs of landmarks,
                          pairedLM = landmark_pairs,
                          ## Enforcing symmetry using this info
                          reflect = TRUE)

# extract the Principle components and store as a data frame
PCscores <- procrustes_data$PCscore_sym
PCscores <- as.data.frame(PCscores)


## Getting the PC scores (the trait space)
traitspace <- PCscores

## to visualise by dietary category
diet_table <- read_csv("files_species_diet.csv")

## Creating a table
## Using the specimens names from the traitspace
diet_table <- data.frame(row.names = rownames(traitspace),
                       ## Add species names 
                       species_name  = diet_table$species_name,
                       ## Add diets
                       diet = diet_table$diet_broad)

## Make sure your table and the traitspace rownames matches (so that the first specimen in your table is the first in your traitspace, etc...)
diet_table <- diet_table[match(rownames(traitspace), rownames(diet_table)), ]

############################ PLOT OF PHYLOGENY WITH THE DIETARY GROUPS
## Read the tree
dataTree <- read.nexus("phylogeny.nex")


## Choosing only the species of interest
reduced_tree <- drop.tip(dataTree,
                         ## Selecting the tips to drop (everything that's not ("!"" = is not)
                         ## in the tree tip labels and in the table species names ("%in%" = intersection))
                         tip = dataTree$tip.label[!(dataTree$tip.label %in% diet_table[, "species_name"])])

# make the .nex file into .tre for easy plotting
write.tree(reduced_tree, file = "phylogeny.tre", append = FALSE,
           digits = 10, tree.names = FALSE)

tree <- read.tree("phylogeny.tre")

tree_plot <- ggtree(tree) %<+%
  diet_table

tree_by_diet <- tree_plot  + xlim(0, 270) +
  geom_tiplab(offset = 5, size = 10) +
  geom_tippoint(aes(color = diet), size = 9, shape = "square") +
  theme(legend.position = "right") +
  scale_color_manual(values=c("#2196F3", "#FFC107", "#9C27B0")) +
  labs(color = "Diet") +
  theme_classic(base_size = 30) +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave(tree_by_diet, filename = "/Users/harveyjeffrey/project_data/plots/tree_by_diet.svg", 
       width = 16, height = 9)



## Plotting the results by diet
diet_PC <- ggplot(traitspace, aes(x = PC1, y = PC2, colour = diet_table$diet)) +
  geom_point(shape = 16, size = 8, alpha = 0.7) +
  scale_colour_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c("#2196F3", "#FFC107", "#9C27B0"),
  ) +
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +
  labs(x = "PC1", y = "PC2", colour = "Diet") +
  ggtitle("(A) Whole Skull") +
  theme(plot.title = element_text(face = "bold"))

## save the plot as an svg file for scalability
ggsave(diet_PC, filename = "/Users/harveyjeffrey/project_data/plots/diet_pc.svg", 
       width = 12, height = 9)



#### partitioning the data 
# List of landmarks that can be superimposed onto the mandible
landmark_subset <- c(4, 5, 6, 7, 8, 10, 11, 16, 17, 18, 22, 23, 24, 25, 26, 28, 29) 

# subset these from the whole landmark list
data_landmark_subset <- scans_array[landmark_subset, , ]

# GPA AND PCA on the subsetted data
procrustes_subset <- procSym(data_landmark_subset)

# Getting the PC scores (the trait space)
traitspace_subset <- procrustes_subset$PCscore

# convert traitspace_subset to a data frame so it works in ggplot
traitspace_subset <- as.data.frame(traitspace_subset)

## Plotting the results by diet for subsetted
diet_PC_subset <- ggplot(traitspace_subset, aes(x = PC1, y = PC2, colour = diet_table$diet)) +
  geom_point(shape = 16, size = 8, alpha = 0.7) +
  scale_colour_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c("#2196F3", "#FFC107", "#9C27B0"),
  ) +
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +
  labs(x = "PC1", y = "PC2", colour = "Diet") +
  ggtitle("(B) Parts mechanically involved in feeding") +
  theme(plot.title = element_text(face = "bold"))

## save the plot as an svg file for scalability
ggsave(diet_PC_subset, filename = "/Users/harveyjeffrey/project_data/plots/diet_pc_subset.svg", 
       width = 12, height = 9)

##### FOR WHAT I DEEM TO BE NOT IMPORTANT IN FEEDING (IE. EVERYTHING - THE SUBSET)
#### partitioning the data 
# List of landmarks that can be superimposed onto the mandible
landmark_cranium <- (1:30)[!(1:30 %in% landmark_subset)]

# subset these from the whole landmark list
data_landmark_cranium <- scans_array[landmark_cranium, , ]

# GPA AND PCA on the subsetted data
procrustes_cranium <- procSym(data_landmark_cranium)

# Getting the PC scores (the trait space)
traitspace_cranium <- procrustes_cranium$PCscore

# convert traitspace_subset to a data frame so it works in ggplot
traitspace_cranium <- as.data.frame(traitspace_cranium)

## Plotting the results by diet for the rest of cranium
diet_PC_cranium <- ggplot(traitspace_cranium, aes(x = PC1, y = PC2, colour = diet_table$diet)) +
  geom_point(shape = 16, size = 8, alpha = 0.7) +
  scale_colour_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c("#2196F3", "#FFC107", "#9C27B0"),
  ) +
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +
  labs(x = "PC1", y = "PC2", colour = "Diet") + 
  ggtitle("(C) Frontal and Parietal") +
  theme(plot.title = element_text(face = "bold"))

ggsave(diet_PC_cranium, filename = "/Users/harveyjeffrey/project_data/plots/diet_pc_cranium.svg", 
       width = 12, height = 9)

PCplot_all <- (diet_PC + diet_PC_subset)/(diet_PC_cranium + plot_spacer())
ggsave(PCplot_all, filename = "/Users/harveyjeffrey/project_data/plots/PC_all.svg", 
       width = 24, height = 18)
######################################### TESTING FOR CONVERGENCE
## Calculating all the traitspace distances
specimen_morphodistances <- as.matrix(dist(traitspace))
## load up the tree
dataTree <- read.nexus("phylogeny.nex")
# calculate evolutionary distances
phylo_distances <- cophenetic(dataTree)

##### to assign distances to a species
get_species_name <- function(specimen, diet_table) {
  return(diet_table[specimen, "species_name"])
}

## This loop replaces each morphological distance by it's equivalent phylogenetic distance
specimen_phylodistances <- specimen_morphodistances
for(specimen in rownames(specimen_morphodistances)) {
  ## Find the species name of the specimen
  morphodistances_to_specimen <- specimen_phylodistances[specimen, ]
  ## Convert the specimen names into species names
  names(morphodistances_to_specimen) <- sapply(names(morphodistances_to_specimen), 
                                               get_species_name, diet_table = diet_table)
  ## Replace the morpho distances by their corresponding phylogenetic distances
  specimen_phylodistances[specimen, ] <- phylo_distances[diet_table[specimen, "species_name"],
                                                         ][names(morphodistances_to_specimen)]
}
specimen_phylodistances <- abs(specimen_phylodistances)
## Getting the morphological distance corrected for phylogenetic distance as per Arbuckle et al (2014)
convergence_ratios <- abs(specimen_morphodistances/(1-log(specimen_phylodistances + 0.01)))

## statistical tests
# select the dietary groups
carnivore <-  which(diet_table[, "diet"] == "Carnivore")
insectivore <-  which(diet_table[, "diet"] == "Insectivore")
generalist <-  which(diet_table[, "diet"] == "Generalist")

## Selects the groups in the distance matrix
select.group <- function(group, specimen_morphodistances) {
  ## Note the [group, group] structure for selecting both rows and columns
  my_matrix <- specimen_morphodistances[group, group]
  ## And the upper.tri function to select only the upper triangle
  ## (because distance matrices are mirrors)
  return(c(my_matrix[upper.tri(my_matrix)]))
} 
carnivore_ratios <- select.group(carnivore, convergence_ratios)
generalist_ratios <- select.group(generalist, convergence_ratios)
insectivore_ratios <- select.group(insectivore, convergence_ratios)


## Testing the differences
kruskal.test(list(carnivore_ratios, generalist_ratios, insectivore_ratios))

0## Strength of convergence
## Function for calculating the wheatseaf index
wheatsheaf.index <- function(group, convergence_ratios) {
  ## Calculate the average convergence ratio
  mean_convergence <- mean(convergence_ratios[upper.tri(convergence_ratios)])
  ## Calculate the average group convergence
  group_ratios <- convergence_ratios[group, group]
  mean_group <- mean(group_ratios[upper.tri(group_ratios)])
  ## return the statistic
  return(mean_convergence/mean_group)
}
## 
wheatsheaf.index(carnivore, convergence_ratios)
wheatsheaf.index(generalist, convergence_ratios)
wheatsheaf.index(insectivore, convergence_ratios)

# make a plot
distance <- data.frame(
  group = rep(c("Carnivore", "Generalist", "Insectivore"), 
              c(length(carnivore_ratios), length(generalist_ratios), 
                length(insectivore_ratios))),
  ratios = c(carnivore_ratios, generalist_ratios, insectivore_ratios),
  subset = "whole_skull"
)

# Plot densities using ggplot2
plot_distance <- ggplot(distance, aes(x = ratios, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#2196F3", "#FFC107", "#9C27B0")) +
  theme_classic(base_size = 30) +
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +  
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Morphological Distance", y = "Density", fill = "Diet") +
  ggtitle("(A) Whole Skull") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_distance, filename = "/Users/harveyjeffrey/project_data/plots/morpho_distance_density.svg", 
       width = 12, height = 9)

###### REPEAT THIS FOR SUBSET
## Calculating all the traitspace distances
specimen_morphodistances_subset <- as.matrix(dist(traitspace_subset))
## load up the tree
dataTree <- read.nexus("phylogeny.nex")
# calculate evolutionary distances
phylo_distances <- cophenetic(dataTree)

##### to assign distances to a species
get_species_name <- function(specimen, diet_table) {
  return(diet_table[specimen, "species_name"])
}

## This loop replaces each morphological distance by it's equivalent phylogenetic distance
specimen_phylodistances_subset <- specimen_morphodistances_subset
for(specimen in rownames(specimen_morphodistances_subset)) {
  ## Find the species name of the specimen
  morphodistances_to_specimen <- specimen_phylodistances_subset[specimen, ]
  ## Convert the specimen names into species names
  names(morphodistances_to_specimen) <- sapply(names(morphodistances_to_specimen), 
                                               get_species_name, diet_table = diet_table)
  ## Replace the morpho distances by their corresponding phylogenetic distances
  specimen_phylodistances_subset[specimen, ] <- phylo_distances[diet_table[specimen, "species_name"],
  ][names(morphodistances_to_specimen)]
}
specimen_phylodistances_subset <- abs(specimen_phylodistances_subset)
## Getting the morphological distance corrected for phylogenetic distance as per Arbuckle et al (2014)
convergence_ratios_subset <- abs(specimen_morphodistances_subset/(1-log(specimen_phylodistances_subset + 0.01)))

## statistical tests
# select the dietary groups
carnivore <-  which(diet_table[, "diet"] == "Carnivore")
insectivore <-  which(diet_table[, "diet"] == "Insectivore")
generalist <-  which(diet_table[, "diet"] == "Generalist")

## Selects the groups in the distance matrix
select.group <- function(group, specimen_morphodistances_subset) {
  ## Note the [group, group] structure for selecting both rows and columns
  my_matrix <- specimen_morphodistances_subset[group, group]
  ## And the upper.tri function to select only the upper triangle
  ## (because distance matrices are mirrors)
  return(c(my_matrix[upper.tri(my_matrix)]))
} 
carnivore_ratios_subset <- select.group(carnivore, convergence_ratios_subset)
generalist_ratios_subset <- select.group(generalist, convergence_ratios_subset)
insectivore_ratios_subset <- select.group(insectivore, convergence_ratios_subset)

## Testing the differences
kruskal.test(list(carnivore_ratios_subset, 
                  generalist_ratios_subset, 
                  insectivore_ratios_subset))

## Strength of convergence through the wheatsheaf index
## Function for calculating the wheatseaf index
wheatsheaf.index <- function(group, convergence_ratios_subset) {
  ## Calculate the average convergence ratio
  mean_convergence_subset <- mean(convergence_ratios_subset[upper.tri(convergence_ratios_subset)])
  ## Calculate the average group convergence
  group_ratios_subset <- convergence_ratios_subset[group, group]
  mean_group_subset <- mean(group_ratios_subset[upper.tri(group_ratios_subset)])
  ## return the statistic
  return(mean_convergence_subset/mean_group_subset)
}
## 
wheatsheaf.index(carnivore, convergence_ratios_subset)
wheatsheaf.index(generalist, convergence_ratios_subset)
wheatsheaf.index(insectivore, convergence_ratios_subset)

# make a plot
distance_subset <- data.frame(
  group = rep(c("Carnivore", "Generalist", "Insectivore"), 
              c(length(carnivore_ratios_subset), length(generalist_ratios_subset), 
                length(insectivore_ratios_subset))),
  ratios = c(carnivore_ratios_subset, generalist_ratios_subset, 
             insectivore_ratios_subset),
  subset = "feeding"
)

# Plot densities using ggplot2
plot_distance_subset <- ggplot(distance_subset, aes(x = ratios, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#2196F3", "#FFC107", "#9C27B0")) +
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
                                        # make the background transparent
                                        plot.background = element_blank(),
                                        panel.background = element_blank()) +  
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Morphological Distance", y = "Density", fill = "Diet") +
  ggtitle("(B) Parts mechanically involved in feeding") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_distance_subset, filename = "/Users/harveyjeffrey/project_data/plots/morpho_distance_subset.svg", 
       width = 12, height = 9)

### Repeat for rest of skull
## Calculating all the traitspace distances
specimen_morphodistances_cranium <- as.matrix(dist(traitspace_cranium))
## load up the tree
dataTree <- read.nexus("phylogeny.nex")
# calculate evolutionary distances
phylo_distances <- cophenetic(dataTree)

##### to assign distances to a species
get_species_name <- function(specimen, diet_table) {
  return(diet_table[specimen, "species_name"])
}

## This loop replaces each morphological distance by it's equivalent phylogenetic distance
specimen_phylodistances_cranium <- specimen_morphodistances_cranium
for(specimen in rownames(specimen_morphodistances_cranium)) {
  ## Find the species name of the specimen
  morphodistances_to_specimen <- specimen_phylodistances_cranium[specimen, ]
  ## Convert the specimen names into species names
  names(morphodistances_to_specimen) <- sapply(names(morphodistances_to_specimen), 
                                               get_species_name, diet_table = diet_table)
  ## Replace the morpho distances by their corresponding phylogenetic distances
  specimen_phylodistances_cranium[specimen, ] <- phylo_distances[diet_table[specimen, "species_name"],
  ][names(morphodistances_to_specimen)]
}
specimen_phylodistances_cranium <- abs(specimen_phylodistances_cranium)
## Getting the morphological distance corrected for phylogenetic distance as per Arbuckle et al (2014)
convergence_ratios_cranium <- abs(specimen_morphodistances_cranium/(1-log(specimen_phylodistances_cranium + 0.01)))

## statistical tests
# select the dietary groups
carnivore <-  which(diet_table[, "diet"] == "Carnivore")
insectivore <-  which(diet_table[, "diet"] == "Insectivore")
generalist <-  which(diet_table[, "diet"] == "Generalist")

## Selects the groups in the distance matrix
select.group <- function(group, specimen_morphodistances_cranium) {
  ## Note the [group, group] structure for selecting both rows and columns
  my_matrix <- specimen_morphodistances_cranium[group, group]
  ## And the upper.tri function to select only the upper triangle
  ## (because distance matrices are mirrors)
  return(c(my_matrix[upper.tri(my_matrix)]))
} 
carnivore_ratios_cranium <- select.group(carnivore, convergence_ratios_cranium)
generalist_ratios_cranium <- select.group(generalist, convergence_ratios_cranium)
insectivore_ratios_cranium <- select.group(insectivore, convergence_ratios_cranium)

## Testing the differences
kruskal.test(list(carnivore_ratios_cranium, 
                  generalist_ratios_cranium, 
                  insectivore_ratios_cranium))

## Strength of convergence through the wheatsheaf index
## Function for calculating the wheatseaf index
wheatsheaf.index <- function(group, convergence_ratios_cranium) {
  ## Calculate the average convergence ratio
  mean_convergence_cranium <- mean(convergence_ratios_cranium[upper.tri(convergence_ratios_cranium)])
  ## Calculate the average group convergence
  group_ratios_cranium <- convergence_ratios_cranium[group, group]
  mean_group_cranium <- mean(group_ratios_cranium[upper.tri(group_ratios_cranium)])
  ## return the statistic
  return(mean_convergence_cranium/mean_group_cranium)
}
## 
wheatsheaf.index(carnivore, convergence_ratios_cranium)
wheatsheaf.index(generalist, convergence_ratios_cranium)
wheatsheaf.index(insectivore, convergence_ratios_cranium)

# make a plot
distance_cranium <- data.frame(
  group = rep(c("Carnivore", "Generalist", "Insectivore"), 
              c(length(carnivore_ratios_cranium), length(generalist_ratios_cranium), 
                length(insectivore_ratios_cranium))),
  ratios = c(carnivore_ratios_cranium, generalist_ratios_cranium, 
             insectivore_ratios_cranium),
  subset = "basicranium"
)

# Plot densities using ggplot2
plot_distance_cranium <- ggplot(distance_cranium, aes(x = ratios, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#2196F3", "#FFC107", "#9C27B0")) +
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap between bars and x-axis#
  labs(x = "Morphological Distance", y = "Density", fill = "Diet") +
  ggtitle("(C) Frontal and Parietal") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_distance_cranium, filename = "/Users/harveyjeffrey/project_data/plots/morpho_distance_cranium.svg", 
       width = 12, height = 9)

distance_all <- (plot_distance + plot_distance_subset)/(plot_distance_cranium + plot_spacer())
ggsave(distance_all, filename = "/Users/harveyjeffrey/project_data/plots/distance_all.svg",
       width = 24, height = 18)

### CONVERGENCE FOR EVERYTHING
all_distances <- rbind(distance, distance_cranium, distance_subset)

anova_distance <- lm(ratios ~ subset, data = all_distances)
anova(anova_distance)

## TUKEY HSD
HSD.test(anova_distance, "subset", console=TRUE)

plot_all_distance <- ggplot(all_distances, aes(x = ratios, fill = subset)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#009688", "#FFDB58", "#3F51B5"),
    limits = c("basicranium", "feeding", "whole_skull"),
    labels = c("Frontal and Parietal", "Feeding parts", "Whole Skull")
  )+
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove gap between bars and x-axis#
  labs(x = "Morphological Distance", y = "Density", fill = "Subset") 

ggsave(plot_all_distance, filename = "/Users/harveyjeffrey/project_data/plots/distance_by_subset.svg",
       width = 12, height = 9)
######################################### CENTROID/SHAPE ANALYSIS

# format the tree so it works in MCMCglmm
dataTreeNode <- makeNodeLabel(dataTree, method = "number") # Add node labels 
INtree <- inverseA(dataTreeNode, nodes="TIPS") #Make matrix of your labelled tree 
Ainv <- INtree$Ainv # Extract the matrix for use in MCMCglmm
prior1 <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002)))


# whole skull
centroid_sizes <- procrustes_data$size
centroid_df <- data.frame(subset = "whole",
                          centroid_size = centroid_sizes,
                             diet = select(diet_table, diet),
                          species = select(diet_table, species_name))
centroid_mcmc <- MCMCglmm(centroid_size ~ diet,
                          random = ~species_name,
                          data = centroid_df,
                          nitt = 100000,
                          verbose = TRUE,
                          prior = prior1, 
                          ginverse = list(species_name = Ainv))


# just partitioned feeding bits
centroid_sizes_subset <- procrustes_subset$size
centroid_subset_df <- data.frame(subset = "jaw",
                                 centroid_size = centroid_sizes_subset,
                                    diet = select(diet_table, diet),
                                 species = select(diet_table, species_name))
centroid_subset_mcmc <- MCMCglmm(centroid_size ~ diet,
                          random = ~species_name,
                          data = centroid_subset_df,
                          nitt = 100000,
                          verbose = TRUE,
                          prior = prior1, 
                          ginverse = list(species_name = Ainv))

# rest of cranium
centroid_sizes_cranium <- procrustes_cranium$size
centroid_cranium_df <- data.frame(subset = "basocranium", 
                                  centroid_size = centroid_sizes_cranium,
                                 diet = select(diet_table, diet),
                                 species = select(diet_table, species_name))
centroid_cranium_mcmc <- MCMCglmm(centroid_size ~ diet,
                                 random = ~species_name,
                                 data = centroid_cranium_df,
                                 nitt = 100000,
                                 verbose = TRUE,
                                 prior = prior1, 
                                 ginverse = list(species_name = Ainv))

summary(centroid_mcmc)
summary(centroid_subset_mcmc)
summary(centroid_cranium_mcmc)


### Centroid graphs
# Plot densities using ggplot2
plot_centroid <- ggplot(centroid_df, aes(x = diet, y = centroid_size, fill = diet)) +
  geom_boxplot(position = "dodge", color = "black", size = 0.5, width = 0.8,
               show.legend = FALSE)  + 
  labs(x = "Diet", y = "Centroid Size") +
  scale_fill_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c( "#2196F3", "#FFC107", "#9C27B0"),
  ) +
  scale_y_continuous(expand = c(0, 0))  + # Remove gap between bars and x-axis
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank(),
        ) +  
  ggtitle("(A) Whole Skull") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_centroid, filename = "/Users/harveyjeffrey/project_data/plots/centroid_whole.svg", 
       width = 12, height = 9)


plot_centroid_subset <- ggplot(centroid_subset_df, aes(x = diet, y = centroid_size, fill = diet)) +
  geom_boxplot(position = "dodge", color = "black", size = 0.5, width = 0.8,
               show.legend = FALSE)  + 
  labs(x = "Diet", y = "Centroid Size") +
  scale_fill_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c( "#2196F3", "#FFC107", "#9C27B0"),
  ) +
  scale_y_continuous(expand = c(0, 0))  + # Remove gap between bars and x-axis
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +  
  ggtitle("(B) Parts mechanically involved in feeding") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_centroid_subset, filename = "/Users/harveyjeffrey/project_data/plots/centroid_subset.svg", 
       width = 12, height = 9)

# basocranium
plot_centroid_cranium <- ggplot(centroid_cranium_df, aes(x = diet, y = centroid_size, fill = diet)) +
  geom_boxplot(position = "dodge", color = "black", size = 0.5, width = 0.8,
               show.legend = FALSE)  + 
  labs(x = "Diet", y = "Centroid Size") +
  scale_fill_manual(
    limits = c("Carnivore", "Generalist", "Insectivore"),
    values = c( "#2196F3", "#FFC107", "#9C27B0"),
  ) +
  scale_y_continuous(expand = c(0, 0))  + # Remove gap between bars and x-axis
  theme_classic(base_size = 30) + 
  theme(axis.title = element_text(face="bold"), 
        # make the background transparent
        plot.background = element_blank(),
        panel.background = element_blank()) +  
  ggtitle("(C) Frontal and Parietal") +
  theme(plot.title = element_text(face = "bold"))

ggsave(plot_centroid_cranium, filename = "/Users/harveyjeffrey/project_data/plots/centroid_cranium.svg", 
       width = 12, height = 9)

centroid_all <- (plot_centroid + plot_centroid_subset)/(plot_centroid_cranium + plot_spacer())

ggsave(centroid_all, filename= "/Users/harveyjeffrey/project_data/plots/centroid_all.svg",
       width = 24, height = 18)
############################################### DISPARITY ANALYSIS
# make the diet table into a matrix that is readable by dispRity
diet_matrix <- as.matrix(diet_table["diet"])

## create a disparity object
traitspace_disparity <- custom.subsets(traitspace, group = diet_matrix)

## bootstrap this
traitspace_disparity_bs <- boot.matrix(traitspace_disparity)

## Measuring disparity
disparity_groups <- dispRity(traitspace_disparity_bs, metric = c(sum, variances))

## Summarizing the results
summary(disparity_groups)

# disparity_groups <- get.disparity(disparity_groups)
# disparity_groups <- as.data.frame(disparity_groups)

plot(disparity_groups)

###### FOR SUBSET
traitspace_subset_disparity <- custom.subsets(traitspace_subset, group = diet_matrix)

## bootstrap this
traitspace_subset_disparity_bs <- boot.matrix(traitspace_subset_disparity)

## Measuring disparity
disparity_groups_subset <- dispRity(traitspace_subset_disparity_bs, metric = c(sum, variances))

## Summarizing the results
summary(disparity_groups_subset)


##### FOR EVERYTHING -SUBSET
traitspace_cranium_disparity <- custom.subsets(traitspace_cranium, group = diet_matrix)

## bootstrap this
traitspace_cranium_disparity_bs <- boot.matrix(traitspace_cranium_disparity)

## Measuring disparity
disparity_groups_cranium <- dispRity(traitspace_cranium_disparity_bs, metric = c(sum, variances))

## Summarizing the results
summary(disparity_groups_cranium)

######## Make a plot of disparity scores
 disparity_data <- read_csv("disparity_data.csv")

disparity_plot <- ggplot(disparity_data, aes(x = subset, fill = diet)) +
  geom_boxplot(
    stat = "identity",
    aes(lower = lower_quartile,
        upper = upper_quartile,
        middle = median,
        ymin   = min, # optional
        ymax   = max)) +
  labs(x = "", y = "Morphological Disparity", fill = "Diet") +  
  scale_fill_manual(
    values = c("#2196F3", "#FFC107", "#9C27B0"),
    labels = c("Carnivore", "Generalist", "Insectivore")) +
  scale_x_discrete(labels = c("Whole Skull", "Frontal and\nParietal", "Feeding Parts")) +
  theme_classic(base_size = 30) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"), 
              # make the background transparent
              plot.background = element_blank(),
              panel.background = element_blank())
ggsave(disparity_plot, filename = "../plots/disparity_plot.svg",
       width = 12, height = 9)


### statsitcal test for disparity
disparity_model <- lm(median ~ subset + diet,
                      data = disparity_data)

anova(disparity_model)