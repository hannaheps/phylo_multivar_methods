library("ape") # dealing with phylogenies
library("phytools") # phylogengenetic methods
library("OUwie") # modeling trait evolution
library("castor") # generating a random transition matrix

# Load HR2015 phylogeny
coral_phy <- ape::read.tree("../input/huang_roy_molecular_r2.newick")

# -- Simulate Host Traits --

# Generate random equal-rates Q matrix with two states
r.Q <- castor::get_random_mk_transition_matrix(
    Nstates = 2,
    rate_model = "ER"
)

# Simulate random discrete HOST trait across tree (root. value)
s.trait_history <- phytools::sim.history(tree = coral_phy, Q = r.Q, anc = "1")

# Simulate uncorrelated continuous HOST trait (evolving via BM)
host_uncorrelated_cont <- OUwie::OUwie.sim(
    phy = s.trait_history,
    simmap.tree = TRUE,
    alpha = c(1e-10, 1e-10),
    sigma.sq = c(1, 1),
    theta0 = 0, theta = c(0, 0)
)


host_df <- data.frame(disc_trait_ER = s.trait_history$states, host_cont_trait = host_uncorrelated_cont$X)

# -- Simulate Microbiome --

# Simulate microbiome

n_microbes <- 1000

raw_microbiome_df <- data.frame(Species = coral_phy$tip.label)

for (i in 1:n_microbes) {
    tmp_microbe <- OUwie::OUwie.sim(phy = s.trait_history, simmap.tree = TRUE, alpha = c(1e-10, 1e-10), sigma.sq = c(1, 1), theta0 = 0, theta = c(0, 0))
    tmp_colname <- paste("microbe", i, sep = "_")
    raw_microbiome_df[[tmp_colname]] <- tmp_microbe$X
    if (i == n_microbes) {
        system("terminal-notifier -title Rstudio -message 'Simulated microbiomes have been generated'")
    }
}


# -- Output --

# Export host data
host_df <- cbind(Species.name = rownames(host_df), host_df)
write.csv(x = host_df, file = "../output/simulated_host_traits.csv", row.names = FALSE)

# Export microbiome data
write.csv(
  x = raw_microbiome_df,
  file = paste("../output/raw_sim_microbiome_",n_microbes,"microbes.csv", sep=""),
  row.names = FALSE
)



# # Combine host and microbe trait data
# 
# library("tibble")
# raw_microbiome_df <- tibble::column_to_rownames(.data = raw_microbiome_df, var = "Species")
# host_and_mbiome_df <- merge(host_df, raw_microbiome_df, by = "row.names")
# host_and_mbiome_df <- tibble::column_to_rownames(.data = host_and_mbiome_df, var = "Row.names")
# unloadNamespace("tibble")
# 
# write.csv(x = host_and_mbiome_df, file = "../output/simulated_mbiome_and_host.csv")
