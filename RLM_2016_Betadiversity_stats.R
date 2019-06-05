
#################################################################################
# This script is the beta diversity analysis for Poc_2016:
# Here, I calculate between and within group distances and 
# statistics for PERMANOVA and PERMDISP tests.
#
# Created by Rebecca Maher
# Created on 10/2/18
# Edited on 03/04/19
#################################################################################

## clear workspace------------------------
rm(list=ls())

# load libraries
library('phyloseq'); packageVersion('phyloseq')
library('vegan'); packageVersion('vegan')

# set working directory-------------------
#setwd("~/Box Sync/RAPID-analysis/")

## functions----------------------
pairwise.adonis.dm <- function(x,factors,stratum=NULL,p.adjust.m="bonferroni",perm=999){
  
  library(vegan)
  
  if(class(x) != "dist") stop("x must be a dissimilarity matrix (dist object)")
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    sub_inds <- factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))
    
    resp <- as.matrix(x)[sub_inds,sub_inds]
    
    ad = adonis(as.dist(resp) ~
                  
                  factors[sub_inds], strata=stratum[sub_inds], permutations=perm);
    
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
}

## Data Analysis----------------------------

# load the rarefied OTU table with mapping file with physiological data
qd <- load("~/data/RLM_2016_phyloseq_object.RData")

# Log-transform OTU-table
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

# Calculate all four distance matrices: weighted unifrac, unweighted unifrac, bray curtis, binary jaccard
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bc <- phyloseq::distance(qd, method = "bray")
qd_bj <- distance(qd, method = "jaccard", binary =TRUE)

ord_wu <- ordinate(qd, "PCoA", distance = qd_wu)
ord_un <- ordinate(qd, "PCoA", distance = qd_un)
ord_bc <- ordinate(qd, "PCoA", distance = qd_bc)
ord_bj <- ordinate(qd, "PCoA", distance = qd_bj)

# PERMANOVA's with Adonis -  Permutational Multivariate Analasis of Variance
# Make a data frame from the sample_data
sampledf <- data.frame(sample_data(qd))

# Adonis individual tests
adonis(qd_bc ~ temp, data = sampledf)
adonis(qd_bj ~ temp, data = sampledf)
adonis(qd_wu ~ temp, data = sampledf)
adonis(qd_un ~ temp, data = sampledf)

adonis(qd_bc ~ nutrient, data = sampledf)
adonis(qd_bj ~ nutrient, data = sampledf)
adonis(qd_wu ~ nutrient, data = sampledf)
adonis(qd_un ~ nutrient, data = sampledf)

adonis(qd_bc ~ corallivory, data = sampledf)
adonis(qd_bj ~ corallivory, data = sampledf)
adonis(qd_wu ~ corallivory, data = sampledf)
adonis(qd_un ~ corallivory, data = sampledf)

# Adonis test for between group diversity, with full formula
adonis(qd_bc ~ temp*corallivory*nutrient, data = sampledf)
adonis(qd_bj ~ temp*corallivory*nutrient, data = sampledf) 
adonis(qd_wu ~ temp*corallivory*nutrient, data = sampledf) 
adonis(qd_un ~ temp*corallivory*nutrient, data = sampledf) 


# PERMDISP with betadisper - Multivariate Homogeneity of group dispersions
anova(betadisper(qd_bc, sampledf$corallivory, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$corallivory, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$corallivory, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$corallivory, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$nutrient, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$nutrient, bias.adjust = TRUE))

anova(betadisper(qd_bc, sampledf$temp, bias.adjust = TRUE))
anova(betadisper(qd_bj, sampledf$temp, bias.adjust = TRUE))
anova(betadisper(qd_wu, sampledf$temp, bias.adjust = TRUE))
anova(betadisper(qd_un, sampledf$temp, bias.adjust = TRUE))

# Pairwise adonis
pairwise.adonis.dm(qd_bc, sample_data(qd)$nutrient, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_bc, sample_data(qd)$temp, p.adjust.m = "fdr")
pairwise.adonis.dm(qd_bc, sample_data(qd)$corallivory, p.adjust.m = "fdr")
# Pairwise betadisper with fdr correction
p.adjust(permutest(betadisper(qd_un, sampledf$nutrient, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_un, sampledf$corallivory, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')
p.adjust(permutest(betadisper(qd_un, sampledf$temp, 
                              bias.adjust = TRUE), 
                   pairwise=TRUE)$pairwise$permuted, method = 'fdr')

###################################################################################
## Extracting data into .csv format for plotting
#### Between group distances
## Add between group distances into a mapping file for plotting
wu_means <- data.frame(rowMeans(as.matrix(qd_wu)))
un_means <- data.frame(rowMeans(as.matrix(qd_un)))
bc_means <- data.frame(rowMeans(as.matrix(qd_bc)))
bj_means <- data.frame(rowMeans(as.matrix(qd_bj)))
# Make a dataframe
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(wu_means) 
# makes a new dataframe with the specified columnes
raw_distances <- data.frame(SampleID, wu_means, un_means, bc_means,bj_means) 
# Makes metadata into a df to work with
s <- data.frame(sample_data(qd)) 
# Change first column title to SampleID to match distances dataframe
colnames(s)[1] <- "SampleID" 
# merges metadata df and distances df
distances <- merge(raw_distances, s, by = "SampleID") 
colnames(distances)[2:5] <- c("distwu", "distun", "distbc","distbj")
betdistmax <- betdist
write.csv(betdistmax, file = "~/RLM_2016_Betadiversity_between_group_distances.csv")

### Within group distances
# Extract within group distances and export as .csv
# Gives the distances of each group to its centroid by interaction
distbc <- betabc$distances
distbj <- betabj$distances
distun <- betaun$distances
distwu <- betawu$distances

# IMPORTANT STEP: Merge the distance to centroid for each distance measure by SampleID
# Put them into the map file for use in the graphs

# Makes a dataframe of the distance to centroid of one distance measure
# Just necessary to get the SampleID labels
df_bc <- data.frame(betabc$distances) 
# Extracts the SampleIDs from the dataframe
SampleID <- row.names(df_bc) 
# makes a new dataframe with the specified columnes
beta_distances <- data.frame(SampleID, distbc, distbj, distun, distwu) 
# Makes metadata into a df to work with
s <- data.frame(sample_data(qd)) 
# Change first column title to SampleID to match distances dataframe
colnames(s)[1] <- "SampleID" 
# merges metadata df and distances df
withdistmax <- merge(beta_distances, s, by = "SampleID") 
write.csv(withdistmax, file = "~/RLM_2016_Betadiversity_within_group_distances.csv")
