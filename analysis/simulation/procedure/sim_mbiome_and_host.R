set.seed(1989)
# TODO -- Roxygen docstrings for different functions | 'roxygen2'
# TODO -- add table dimension checks AND other test code

library("ape") # dealing with phylogenies
library("phytools") # phylogengenetic methods
library("OUwie") # modeling trait evolution
library("castor") # generating a random transition matrix
library("optparse") # parse arguments python-style

cat("Done loading packages. \n\n")

# command line argument parsing
option_list <- list (
  optparse::make_option(
    c("--sim_number"),
    type = "integer",
    help = "Simulation number (INTEGER) will produce S##_[sim_info].zip"
  ),
  optparse::make_option(
    c("--n_ASV"),
    type = "integer",
    action = "store",
    help = "Number of microbes to simulate (INTEGER)"
  ),
  optparse::make_option(
    c("--microbe_sim_type"),
    type = "character",
    help = "Feature simulation type (CHARACTER). OPTIONS: random, bm, ou"
  ),
  optparse::make_option(
    c("--n_simulations"),
    type = "integer",
    help = "Number of tables to produce (INTEGER)"
  ),
  optparse::make_option(
    c("--out_dir"),
    type = "character",
    help = "Location where folder should be placed, example: '../output/' (CHARACTER)"
  ),
  optparse::make_option(
    c("--effect_size"),
    type = "double",
    help = "Effect size when OU model is selected (DOUBLE)"
  ),
  optparse::make_option(
    c("--host_trait_type"),
    type = "character",
    # default = FALSE,
    help = "Generate continuous host trait folder (cont) or discrete host trait (disc) (CHARACTER)"
  ),
  optparse::make_option(
    c("--feature_cont_trait_prop_corr"),
    type = "double",
    default = 1,
    help = "Set proportion of ASVs correlated -- should be between 0 and 1 (DEFAULT/1) (DOUBLE)"
  ),
  optparse::make_option(
    c("--feature_cont_trait_corr"),
    type = "double",
    help = "Set mean covariance for continuous host trait with feature traits (DOUBLE)"
  ),
  optparse::make_option(
    c("--feature_cont_trait_corr_sd"),
    type = "double",
    help = "Set standard deviation of covariance for cont. host trait (DOUBLE)"
  ),
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)


if(opt$microbe_sim_type == "ou"){
  if(is.null(opt$effect_size)){
    stop("ERROR: Need to specify effect size if OU microbe model selected")
  } else {}
} else {}

if(opt$host_trait_type == "cont"){
  if (opt$feature_cont_trait_prop_corr < 0 || opt$feature_cont_trait_prop_corr > 1){
    stop("ERROR: Proportion correlated should between 0 and 1")
  } else {}
  if(is.null(opt$feature_cont_trait_corr) || is.null(opt$feature_cont_trait_corr_sd)){
    stop("ERROR: Need to specify correlation value AND stdev")
  } else {}
} else {}

# check if all options are specified, if not stop and help
if(any(c(is.null(opt$n_ASV), is.null(opt$n_simulations), is.null(opt$microbe_sim_type), is.null(opt$sim_number)), is.null(opt$out_dir)) == TRUE) {
  cat("ERROR: One or more options haven't been specified.\n")
  cat(paste("Simulation ID missing?", is.null(opt$sim_number)))
  cat(paste("# of ASVs missing?", is.null(opt$n_ASV)), "\n")
  cat(paste("# of sims missing?", is.null(opt$n_simulations)), "\n")
  cat(paste("Simulation type missing?", is.null(opt$microbe_sim_type)), "\n")
  cat(paste("Output dir missing?", is.null(opt$out_dir)), "\n\n")
  optparse::print_help(opt_parser)
  stop("Please specify fill all options (see error note above)")
}

sim_number <- opt$sim_number
n_ASV <- opt$n_ASV
n_simulations <- opt$n_simulations
microbe_sim_type <- opt$microbe_sim_type
out_dir <- opt$out_dir
effect_size <- opt$effect_size
host_trait_type <- opt$host_trait_type
cont_trait_corr <- opt$feature_cont_trait_corr
cont_trait_corr_sd <- opt$feature_cont_trait_corr_sd
cont_trait_prop_corr <- opt$feature_cont_trait_prop_corr

# components of information for archive file
simulation_id = paste("S", sprintf('%02d', sim_number), sep="")
host_info = "host_ER_BM"
microbe_info = paste("microbe", microbe_sim_type, n_ASV, sep = "_")
tree_info = "coral_tree"
simulation_info = paste("x", n_simulations, sep="")

archive_name = paste(simulation_id, host_info, microbe_info, tree_info, simulation_info, sep = "_")
output_zip = paste(out_dir, archive_name, ".zip", sep="")


cat("\nGenerating a zip folder of simulated microbial data files with these attributes:\n")
cat(paste("Simulation ID:", opt$n_ASV, "\n"))
cat(paste("Total ASVs:", opt$n_ASV, "\n"))
cat(paste("Total tables:", opt$n_simulations, "\n"))
cat(paste("Continuous Similuation?"))
if(host_trait_type == "cont") {
  cat(paste("Correlation with host:", opt$))
}
cat(paste("Simulation type:", opt$microbe_sim_type, "\n"))
if(exists("opt$effect_size")){
  cat(paste("Effect size:", opt$effect_size,"\n"))
}
cat(paste("Output archive:", output_zip), "\n--\n\n")



archive_base_name <- gsub(".*/([^.]+)\\.zip$", "\\1", output_zip) # for zip


tmp_output_folder <- paste("./", archive_base_name, "/", sep = "")

unlink(tmp_output_folder, recursive = TRUE) # delete temporary output folder (in case of prior failed run)

dir.create(tmp_output_folder) # create tmp output folder

# Define functions for each simulation type
random_uncorrelated <- function(n_ASV, tmp_microbes) {
  for (j in 1:n_ASV) {
      tmp_trait <- rnorm(n = n_species, mean = 0, sd = 10)
      tmp_colname <- paste("microbe", j, sep = "_")
      tmp_microbes[[tmp_colname]] <- tmp_trait
  }
  return(tmp_microbes)
}

bm_uncorrelated <- function(n_ASV, tmp_microbes) {
  for (j in 1:n_ASV) {
      tmp_trait <- OUwie::OUwie.sim(phy = s.trait_history, simmap.tree = TRUE, alpha = c(1e-10, 1e-10), sigma.sq = c(1, 1), theta0 = 0, theta = c(0, 0))
      tmp_colname <- paste("microbe", j, sep = "_")
      tmp_microbes[[tmp_colname]] <- tmp_trait$X
  }
  return(tmp_microbes)
}

ou_correlated <- function(n_ASV, tmp_microbes) {
  # prop_affected <- 0.12 -- WE ARE MAKING 
  # effect_size # for transforming theta
  # microbe_effect_df <- data.frame(matrix(nrow = 0, ncol = 2))
  # colnames(microbe_effect_df) <- c("ASV_id", "effect_size")
  for (j in 1:n_ASV) {
      tmp_colname <- paste("microbe", j, sep = "_")
      # roll <- runif(n = 1, min = 0, max = 1) # roll for a random number between 0 and 1
      # if (roll <= prop_affected) {
      #   tmp_alpha <- 1e-2
      #   microbe_effect_df[nrow(microbe_effect_df) + 1,] <- c(tmp_colname, effect_size)
      # } else if (roll > prop_affected) {
      #   tmp_alpha <- 1e-10
      #   microbe_effect_df[nrow(microbe_effect_df) + 1,] <- c(tmp_colname, 0)
      # }
      tmp_alpha <- 1e-2
      # microbe_effect_df[nrow(microbe_effect_df) + 1,] <- c(tmp_colname, effect_size)
      theta_a <- rnorm(n = 1, mean = 0, sd = 1)
      theta_b <- theta_a*effect_size
      thetas <- sample(c(theta_a, theta_b))
  
      tmp_trait <- OUwie::OUwie.sim(phy = s.trait_history, simmap.tree = TRUE, alpha = c(tmp_alpha, tmp_alpha), sigma.sq = c(1, 1), theta0 = 0, theta = thetas)
      tmp_microbes[[tmp_colname]] <- tmp_trait$X
  }
  return(tmp_microbes)
  # assign("microbe_effect_df", microbe_effect_df, envir = .GlobalEnv)
}

generate_cov_matrix <- function(n_ASVs, host_ASV_cov, host_ASV_cov_sd, prop_corr) {
  # TODO - parameters
  #   percent correlated
  #   degree of correlation (normal distribution around host_ASV_cov, with some standard deviation)
  #   -- rnorm(host_ASV_cov, host_ASV_cov_sd)

  # Initializer n+1 x n+1 matrix with zeroes, where n+1 is the host and number of ASVs 
  cov_matrix <- matrix(0, n_ASVs+1, n_ASVs+1)
  
  # Create degrees of correlated that are normally distributed around host_ASV_cov with stdev host_ASV_cov_sd
  total_correlated_ASVs <- round(prop_corr * n_ASVs)
  total_uncorrelated_ASVs <- n_ASVs - total_correlated_ASVs
  covariance_vec <- c(rnorm(n = total_correlated_ASVs, mean=host_ASV_cov, sd=host_ASV_cov_sd), rep(0, total_uncorrelated_ASVs))
  
  # Set the covariance of all traits with host to 'host_ASV_cov'
  cov_matrix[1,] <- c(1, covariance_vec)
  cov_matrix[,1] <- c(1, covariance_vec)
  
  # Set variances (diagonals) to 1
  diag(cov_matrix) <- 1
  
  return(cov_matrix)
}

# Load HR2015 phylogeny
coral_phy <- ape::read.tree("../input/huang_roy_molecular_r2.newick")
n_species <- length(coral_phy$tip.label)


cat("Simulating host and microbiome...\n")
pb = txtProgressBar(min = 0, max = n_simulations, initial = 0, style = 3) # create progress bar

# if host_trait_type == FALSE
if (host_trait_type == "disc") {
  for (i in 1:n_simulations){
    sim_iteration = sprintf('%04d', i)
    # Simulate host
    # Generate random equal-rates Q matrix with two states

    n.states <- 2
    transition_model <- "ER"

    r.Q <- castor::get_random_mk_transition_matrix(
      Nstates = n.states,
      rate_model = transition_model
    )
    colnames(r.Q) <- c('A','B')
    rownames(r.Q) <- c('A','B')

    # Simulate random discrete HOST trait across tree (root. value)
    s.trait_history <- phytools::sim.history(tree = coral_phy, Q = r.Q, anc = "A", message = FALSE)

    # Simulate uncorrelated continuous HOST trait (evolving via BM)
    host_uncorrelated_cont <- OUwie::OUwie.sim(
      phy = s.trait_history,
      simmap.tree = TRUE,
      alpha = c(1e-10, 1e-10),
      sigma.sq = c(1, 1),
      theta0 = 0, theta = c(0, 0)
    )

    tmp_host_df <- data.frame(disc_trait_ER = s.trait_history$states, host_cont_trait = host_uncorrelated_cont$X)
    tmp_host_df <- cbind('#SampleID' = rownames(tmp_host_df), tmp_host_df)
    
    # Export host table
    write.table(
      x = tmp_host_df,
      file = paste(tmp_output_folder, simulation_id, "_", sim_iteration, "_raw_host_traits.tsv", sep = ""),
      row.names = FALSE,
      sep="\t"
    )
    
    # Simulate microbiome
    tmp_microbes <- data.frame('#SampleID' = coral_phy$tip.label, check.names = FALSE)
    
    if (microbe_sim_type == "random") {
      random_uncorrelated(n_ASV, tmp_microbes)
    } else if (microbe_sim_type == "bm") {
      bm_uncorrelated(n_ASV, tmp_microbes)
    } else if (microbe_sim_type == "ou") {
      ou_correlated(n_ASV, tmp_microbes)
    } else {
      stop("Invalid microbe_sim_type, should be random, bm, or ou")
    }

    # Export microbiome
    write.table(
      x = tmp_microbes,
      file = paste(tmp_output_folder, simulation_id, "_", sim_iteration, "_raw_microbial_traits.tsv", sep = ""),
      row.names = FALSE,
      sep = "\t"
    )

    # Export microbe effect df to see what microbes should be affected
    # ONLY IF OU
    if (exists("microbe_effect_df")){
      write.table(
        x = microbe_effect_df,
        file = paste(tmp_output_folder, simulation_id, "_", sim_iteration, "_raw_microbial_trait_effect_size_dict.tsv", sep = ""),
        row.names = FALSE,
        sep = "\t"
      )
    }
    
    setTxtProgressBar(pb,i) # add to progress bar
  }
} else if (host_trait_type == "cont"){
  for (i in 1:n_simulations){
    sim_iteration = sprintf('%04d', i)
    
    tmp_vcv <- generate_cov_matrix(n_ASVs = n_ASV, host_ASV_cov = cont_trait_corr, host_ASV_cov_sd =  cont_trait_corr_sd, prop_corr = cont_trait_prop_corr)
    tmp_cov_trait <- data.frame(phytools::sim.corrs(tree = coral_phy, vcv = tmp_vcv))
    vcv_colnames <- c("host_trait")
    for (i in 1:n_ASV){
      vcv_colnames[i+1] <- paste("microbe", i, sep = "_") # i+1 to account for host
    }
    colnames(tmp_cov_trait) <- vcv_colnames
    
    tmp_host_and_microbe_df <- cbind('#SampleID' = coral_phy$tip.label, tmp_cov_trait)
    
    write.table(
      x = tmp_host_and_microbe_df,
      file = paste(tmp_output_folder, simulation_id, "_", sim_iteration, "_host_and_microbe_traits.tsv", sep=""),
      row.names = FALSE,
      sep = "\t"
    )

    write.table(
      x = tmp_vcv,
      file = paste(tmp_output_folder, simulation_id, "_", sim_iteration, "_vcv.tsv", sep=""),
      row.names = FALSE,
      sep = "\t",
      headings # FIXME - I don't remember the syntax
    )

    # TODO -- write table for vcv matrix; rename columns to vcv_colnames
    setTxtProgressBar(pb,i) # add to progress bar
  }
} else {}

close(pb) # end progress bar
# -- Output --

# # Export host data
# host_df <- cbind(Species.name = rownames(host_df), host_df)
# write.table(x = host_df, file = "../output/simulated_host_traits.csv", row.names = FALSE, sep = "\t")

# Export and report microbiome data

cat("Archiving simulated data...\n")
zip(zipfile = output_zip, files = tmp_output_folder, extras = "-j")

cat("\nUncompressed folder size\n")
system(paste('du -h', tmp_output_folder)) # Uncompressed folder size

cat("\nCompressed archive size\n")
system(paste('du -sh', output_zip))

cat("\nDeleting uncompressed folder...\n")
unlink(tmp_output_folder, recursive = TRUE) # delete temporary output folder



# # Combine host and microbe trait data
# 
# library("tibble")
# raw_microbiome_df <- tibble::column_to_rownames(.data = raw_microbiome_df, var = "Species")
# host_and_mbiome_df <- merge(host_df, raw_microbiome_df, by = "row.names")
# host_and_mbiome_df <- tibble::column_to_rownames(.data = host_and_mbiome_df, var = "Row.names")
# unloadNamespace("tibble")
# 
# write.csv(x = host_and_mbiome_df, file = "../output/simulated_mbiome_and_host.csv")
