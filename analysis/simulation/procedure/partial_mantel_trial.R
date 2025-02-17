##Code for Partial Mantel test##
##Multivariate Methods##

#Load libraries
library(vegan)

#Set data objects
tree <- ape::read.tree("../input/huang_roy_molecular_r2.newick")
traits <- read.csv("../output/simulated_host_traits.csv", head = T, row.names = 1)
micro <- read.csv("../output/raw_sim_microbiome_1000microbes.csv", head = T, row.names = 1)

#Build distance matrices
host.phylo.dist <- as.dist(cophenetic(tree))
host.trait.dist <- vegdist(traits, method = "bray")
micro.comm.dist <- vegdist(micro, method = "bray")

mantel.partial(host.phylo.dist, host.trait.dist, micro.comm.dist, method = "pearson", permutations = 999)

##Works! Some warnings on the negative values


###Example with a for loop using quick simulations

# Example: Load or generate the phylogenetic distance matrix
tree <- ape::rtree(10)  # Random tree with 10 taxa
phylo_dist_matrix <- as.dist(cophenetic(tree))

# Example: Create lists of trait and environmental distance matrices
trait_dist_list <- list(
  trait1 = vegdist(matrix(rnorm(100), nrow=10)),  # Random trait distance matrix 1
  trait2 = vegdist(matrix(rnorm(100), nrow=10))   # Random trait distance matrix 2
)

env_dist_list <- list(
  env1 = vegdist(matrix(rnorm(100), nrow=10)),  # Random environmental distance matrix 1
  env2 = vegdist(matrix(rnorm(100), nrow=10))   # Random environmental distance matrix 2
)

# Initialize an empty dataframe to store results
results <- data.frame(
  Trait = character(),
  Environment = character(),
  Mantel_R = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through combinations of trait and environmental matrices
for (trait_name in names(trait_dist_list)) {
  for (env_name in names(env_dist_list)) {
    
    # Extract current trait and environmental distance matrices
    trait_dist <- trait_dist_list[[trait_name]]
    env_dist <- env_dist_list[[env_name]]
    
    # Perform the partial Mantel test
    mantel_result <- mantel.partial(phylo_dist_matrix, trait_dist, env_dist, method = "pearson", permutations = 999)
    
    # Store results
    results <- rbind(results, data.frame(
      Trait = trait_name,
      Environment = env_name,
      Mantel_R = mantel_result$statistic,
      P_value = mantel_result$signif
    ))
  }
}

# View results
print(results)

