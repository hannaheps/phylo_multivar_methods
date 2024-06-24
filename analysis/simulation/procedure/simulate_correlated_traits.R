library(ape)
library(phytools)

#TODO: add command line interface

modeltest_tree <- read.tree('../input/huang_roy_molecular_r2.newick')
ntraits <- 1000
output_dir <- "../output"

#Simulate 2 uncorrelated traits with brownian motion
#modeltest_BM_X <- fastBM(modeltest_tree, nsim=ntraits)
#write.csv(modeltest_BM_X, paste0('../output/modeltestdata_BM_',ntraits,'.csv'), row.names = TRUE)

#Simulate 2 uncorrelated traits with Ornstein-Uhlenbeck dynamics, alpha=1.0
#modeltest_OU_a1_X <- fastBM(modeltest_tree, nsim=ntraits, alpha=1,theta=0)
#write.csv(modeltest_OU_a1_X, '../output/modeltestdata_OU_a1.csv', row.names = TRUE)

#Simulate 2 uncorrelated traits with Ornstein-Uhlenbeck dynamics, alpha=2.0
#modeltest_OU_a0.5_X <- fastBM(modeltest_tree, nsim=ntraits, alpha=0.5,theta=0)
#write.csv(modeltest_OU_a1_X, '../output/modeltestdata_OU_a0.5.csv', row.names = TRUE)

#Simulate correlated traits with random uniform correlations

#TODO:
#Convert microbiome simulation code into a function

vcv<-clusterGeneration::genPositiveDefMat(5,"unifcorrmat")$Sigma
vcv

#Generate trait data using simulated traits
modeltest_correlated_traits <- sim.corrs(modeltest_tree, vcv)

colnames(modeltest_correlated_traits) <- c("Trait1","Trait2","Trait3","Trait4","Trait5")
head(modeltest_correlated_traits)

#modeltest_correlated_traits$taxon <- row.names(modeltest_correlated_traits)

#To avoid a blank rowname column, make a Tidy-format version of the table
#in which the data holds the row names, then just write to file "without" official rownames
#Hat tip: https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames

tidy_table <- data.frame("Taxon"=rownames(modeltest_correlated_traits),modeltest_correlated_traits)
write.table(tidy_table, paste0("../output/raw_sim_microbiome_",ntraits,"_traits.csv"), row.names=FALSE,sep=",")


#Calculate cophenetic phylogenetic VCV matrix using phytools 
phylogenetic_VCV <- vcvPhylo(modeltest_tree,anc.nodes=FALSE)
head(phylogenetic_VCV)

write.csv(phylogenetic_VCV, "../output/raw_sim_microbiome_tree_VCV.csv", row.names=TRUE)
write.csv(vcv, "../output/raw_sim_microbiome_trait_VCV.csv",row.names=TRUE)
