#################################################################################
# This script is used to generate alpha diversity statistics for Poc 2016:
# Here, I calculate statistics for community richness, evenness,
# and phylogenetic diversity.
#
# Created by Rebecca Maher
# Created on 10/2/18
# Edited on 03/4/19
#################################################################################

## clear workspace------------------------
rm(list=ls())

## functions
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}

# load the rarefied OTU table with mapping file with physiological data
qd <- load("~/data/RLM_2016_phyloseq_object.RData")

##
# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(qd)

richness <- matrix(nrow = nsamp)
row.names(richness) <- sample_names(qd)

evenness <- matrix(nrow =nsamp)
row.names(evenness) <- sample_names(qd)

faithPD <- matrix(nrow = nsamp)
row.names(faithPD) <- sample_names(qd)
##


##
# Options for measures = ("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

# Calculate richness
rich <- as.numeric(as.matrix(subset(estimate_richness(qd, measures = "Chao1"), select = c(1))))
richness[ ,] <- rich
colnames(richness) [1] <- "richness"

# Calculate evenness
even <- as.numeric(as.matrix(estimate_richness(qd, measures = "Simpson")))
evenness[ ,] <- even
colnames(evenness) [1] <- "evenness"

# Calculate Faith's PD
faith <- as.numeric(as.matrix(subset(estimate_pd(qd), select = c(1))))  # estimate_pd is a function assigned at the end of the script
faithPD[ ,] <- faith
colnames(faithPD) [1] <- "faithPD"

# Included the subset in "rich" because the Chao1 measurement outputs two measures per sample (Chao1 and se.chao1)
# and we only want Chao1, so we select for the first column
##


##
# Combine our estimates for richness and evenness into one dataframe
alpha <- cbind(richness, evenness[,1][match(rownames(richness), rownames(evenness))],
               faithPD[,1][match(rownames(richness), rownames(faithPD))])
colnames(alpha) [2] <- "evenness"
colnames(alpha) [3] <- "faithPD"

# Add the sample metadata into this dataframe using the merge() command
# I had to change my column name from X.SampleID to SampleID to match alpha
colnames(sample_data(qd))[1] <- "SampleID" 
s <- data.frame(sample_data(qd))
alphadiv <- cbind(alpha, s)
write.csv(alphadiv, file = "~/alphadiv.csv")
##
