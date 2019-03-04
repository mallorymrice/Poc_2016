########################################################
### Importing microbiome data from Qiime to Phyloseq 
### and rarefying table for Poc_2016 analyses
### By Rebecca Maher
### Edited on 03/04/19
########################################################

rm(list=ls())

library(phyloseq)
library(ape)

##
# Importing data and transformations ##
# Import qiime mapping file, biom otu table, and tree
mapfile = "~/map.txt"
map = import_qiime_sample_data(mapfile)
class(map)

# Import tree
tree = read_tree("~/rep_set.tre")

# This biom file is a result of the pick_open_reference_otu command in qiime1,
# In qiime1, I already removed mitochondria and chloroplasts, removed otus with < 100 occurrences
# and samples with < 1000 counts.
biomfile = "~/otu_table_mc2_w_tax_no_pynast_failures_newnames_clean_o100_s1000_noneg_filtagain.biom"
biom = import_biom(biomfile, parseFunction = parse_taxonomy_default)

qd = merge_phyloseq(map,tree,biom)

##
# Various necessary adjustments to the data (Specific to your data)
# Change rank names from Rank1, Rank2, ... to Kingdom, Phylum, etc.
colnames(tax_table(qd)) = c(k="Kingdom", p="Phylum", c="Class", o="Order",f="Family", g="Genus", s="Species")
# check if it worked 
rank_names(qd)
# Make a numeric factor into a categorical
# In this case, i had temperature as 26 or 29, R considers this numerical but I want to consider it categorical
sample_data(qd)$temp=factor(get_variable(qd,"temp"))

# Setting the order of variables for downstream plotting
qd@sam_data[["interaction"]] <- factor(qd@sam_data[["interaction"]], c("Control","High","Scarred",
                                                                       "NO3-","NH4+","High, scarred",
                                                                       "NO3-, scarred","NH4+, scarred",
                                                                       "NO3-, high","NH4+, high",
                                                                       "NO3-, high, scarred",
                                                                       "NH4+, high, scarred"))
qd@sam_data[["stress"]] <- factor(qd@sam_data[["stress"]], c("none", "single","double","triple"))
##


##
# Setting the seed to 999 ensures that the rarefied biom table is reproducible
min_lib <- min(sample_sums(biom))
rarebiom <- rarefy_even_depth(biom,sample.size = min_lib, verbose = TRUE, replace = FALSE, rngseed = 999)

qd = merge_phyloseq(map,tree,rarebiom)
##


##
# Save the Formal class phyloseq to an external file to load in other scripts
save(qd, file = "/Users/Becca/Box Sync/CHOMPIN/Manuscript_mal16/R-info/qd.RData")

##