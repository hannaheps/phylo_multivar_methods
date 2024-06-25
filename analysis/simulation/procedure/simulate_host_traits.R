library('ape')
library('phytools')

# load phylogeny (HR2015)
coral_phy <- ape::read.tree("../input/huang_roy_molecular_r2.newick")

# simulate continuous trait evolving via a Brownian motion process
s.cont_trait <- phytools::fastBM(coral_phy)

# simulate discrete trait evolving using an equal rates model (ER) and
#   an all-rates-different (ARD) model

### Equal rates
s.disc_trait_ER <- ape::rTraitDisc(coral_phy, model = "ER")
#### re-simulate this discrete trait until the trait is not uniform betweeen sp
####   (e.g., all species have trait A)
while (length(unique(unname(s.disc_trait_ER))) == 1) {
	s.disc_trait_ER <- ape::rTraitDisc(coral_phy, model = "ER")
}

### All-rates-different (Transition rates A -> B (0.1), B -> A (0.01))
s.disc_trait_ARD <- ape::rTraitDisc(coral_phy, model = "ARD", rate = c(0.1,0.01))
#### re-simulate this discrete trait until the trait is not uniform betweeen sp
####   (e.g., all species have trait A)
while (length(unique(unname(s.disc_trait_ARD))) == 1) {
	s.disc_trait_ARD <- ape::rTraitDisc(coral_phy, model = "ARD", rate = c(0.1, 0.01))
} 

# create dataframe with continuous and discrete trait
s.trait_df <- data.frame(Species.name = coral_phy$tip.label,
						   cont_trait = s.cont_trait,
						   disc_trait_ER = s.disc_trait_ER,
						   disc_trait_ARD  =s.disc_trait_ARD)

# write trait data to a csv file
write.csv(x = s.trait_df, file = '../output/simulated_host_traits.csv', row.names = FALSE)
