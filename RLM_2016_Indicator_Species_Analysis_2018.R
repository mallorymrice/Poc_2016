#####################################################
## Indicator Species Analysis for POC_2016 ##
## November 20, 2018
## By Rebecca Maher
####################################################

rm(list=ls())

library("phyloseq")
library("data.table")
library("plyr")
library("dplyr")
library("ggplot2")
library("reshape2")
library("indicspecies")
library("ggnetwork")
library("ape")


# load the rarefied OTU table with mapping file with physiological data
qd <- load("~/data/RLM_2016_phyloseq_object.RData")

# Subset phyloseq object to only include samples with a healing rate.
qd <- prune_samples(sample_data(qd)$healing !="NA", qd)

# Make a binary table
scarred <- as.data.frame(t(otu_table(qd)))
scar <- as.data.frame(ifelse(scarred>0,1,0))
mean(rowSums(scar))

# With the full table
binary <- as.data.frame(t(otu_table(qd)))
binary <- as.data.frame(ifelse(binary>0,1,0))

# Making the named vector for the cluster option
v = sample_data(qd)$nutrient
names(v) = rownames(sample_data(qd))
# Run indicator species analysis
vals <- multipatt(scar, cluster = v, func = "r.g", control = how(nperm =999))
summary(vals, indvalcomp = TRUE)
vals1 <- signassoc(scar, U = NULL, cluster = v, mode = 1, control = how(nperm = 999))

isa <- read.csv(file = "/Users/Becca/Documents/isa.csv")

####################################################################################
# making figure dot plot
qd = transform_sample_counts(qd, function(x) x/sum(x)) # turn into relative abundances

### Reextract your relative abundance data as community structure and taxonomic annotation
comm = as.data.frame(as(object = phyloseq::otu_table(qd), Class = "matrix"))
tax = as.data.frame(as(object = phyloseq::tax_table(qd), Class = "matrix"))
meta = as.data.frame(as(object = phyloseq::sample_data(qd), Class = "matrix"))
data = cbind(tax,comm)

### Melt your data
data$otuid <- rownames(data)
data_melt = melt(data[,c(8:36)], id = c("otuid"))  # Melt data by Family and all samples
tax_melt = data[,c(3:5,36)]
tax_melt$otuid <- rownames(tax_melt)
new_data_melt <- merge(data_melt, tax_melt)

meta_t <- meta[,c(1,3)]
colnames(meta_t)[1] <- "variable"
temp <- merge(new_data_melt, meta_t)


isa <- c("321405","552866","1116027","267843","New.CleanUp.ReferenceOTU235","New.CleanUp.ReferenceOTU252",
         "916277","509548","New.CleanUp.ReferenceOTU19","706432","256116")
temp <- temp[temp$otuid %in% isa,]
temp$Family <- gsub("f__", "", temp$Family)
temp[,6][is.na(temp[,6])] = "Class: Alphaproteobacteria"
temp$otuid <- as.factor(temp$otuid)
levels(temp$otuid)
temp$otuid <- factor(temp$otuid, levels = c("321405","552866","1116027","267843","New.CleanUp.ReferenceOTU235","New.CleanUp.ReferenceOTU252",
                                            "916277","509548","New.CleanUp.ReferenceOTU19","706432","256116"))

# plot
colorblind_pallette = c("#999999", "#E69F00", "#56B4E9", "#ffffff", "#009E73", "#660066", "#FFFF00", "#0072B2", "#CC3300", "#CC79A7","#000000")
p <- ggplot(temp,aes(variable,otuid)) + 
  geom_point(aes(size=value, fill = Family), shape = 21) + 
  facet_grid(. ~ nutrient,scales = "free", space = "free", switch = "y") +
  theme_facet() +
  scale_fill_manual(values = colorblind_pallette) +
  guides(fill = guide_legend(override.aes = list(size=4))) +
  labs(fill = "Family", size = "Relative abundance")
p
