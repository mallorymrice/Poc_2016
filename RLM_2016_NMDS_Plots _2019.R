#################################################
## Code for plotting NMDS figures
## By Rebecca Maher
## Created 6/5/2019
###############################################

library("phyloseq")
library("cowplot")
library("ggthemes")

# load the rarefied OTU table with mapping file with physiological data
qd <- load("~/data/RLM_2016_phyloseq_object.RData")

# Log-transform OTU-tabel
otus_log <- as(otu_table(qd), "matrix")
otus_log <- decostand(otus_log, method = "log")
OTU_log <- otu_table(otus_log, taxa_are_rows = TRUE)
otu_table(qd) <- OTU_log

# Prepare data for plotting
meta <- as.data.frame(sample_data(qd))
table <- as.data.frame(otu_table(qd))
head(meta)

control <- rownames(meta[which(meta[,3] == "Control"),])
nitrate <- rownames(meta[which(meta[,3] == "NO3-"),])
ammon <- rownames(meta[which(meta[,3] == "NH4+"),])

control.df <- table[,control]
nitrate.df <- table[,nitrate]
ammon.df <- table[,ammon]

test.df <- cbind(control.df, nitrate.df, ammon.df)
classes <- c(rep("Control", length(control)), rep("Nitrate", length(nitrate)), rep("Ammonium", length(ammon)))
groups <- classes
df <- test.df
dims <- c(1,2)
ellp.kind <- "ehull"

# For Weighted Unifra
qd_wu <- phyloseq::distance(qd, method = "wunifrac")
object <- metaMDS(qd_wu, k=3, binary = FALSE)

## Weighted Unifrac ehull
mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims) 
points(mds.fig, "sites", pch = 19, col = "#000000", select = control) 
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = nitrate) 
points(mds.fig, "sites", pch = 19, col = "#56B4E9", select = ammon) 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#000000"), lwd = 2, show.groups = "Control") 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#E69F00"), lwd = 2, show.groups = "Nitrate") 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#56B4E9"), lwd = 2, show.groups = "Ammonium")
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#000000"), lwd = 2, show.groups = "Control", lty = 3) 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#E69F00"), lwd = 2, show.groups = "Nitrate", lty = 3)  
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#56B4E9"), lwd = 2, show.groups = "Ammonium", lty = 3) 
#levels(classes) <- c("Control","Nitrate", "Ammonium")
legend("bottomright", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("#000000", "#E69F00", "#56B4E9"), y.intersp = 2)

# For Binary Jaccard
object <- metaMDS(t(df), distance = "jaccard", k=2, binary = TRUE)

mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims)
points(mds.fig, "sites", pch = 19, col = "#000000", select = control)
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = nitrate)
points(mds.fig, "sites", pch = 19, col = "#56B4E9", select = ammon)
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#000000"), lwd = 2, show.groups = "Control")
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#E69F00"), lwd = 2, show.groups = "Nitrate")
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#56B4E9"), lwd = 2, show.groups = "Ammonium")
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#000000"), lwd = 2, show.groups = "Control", lty = 3)
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#E69F00"), lwd = 2, show.groups = "Nitrate", lty = 3)
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#56B4E9"), lwd = 2, show.groups = "Ammonium", lty = 3)
levels(classes) <- c("Control","Nitrate", "Ammonium")
legend("bottomleft", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("#000000", "#E69F00", "#56B4E9"), y.intersp = 2)
