#################################################
## Code for plotting NMDS figures
## By Rebecca Maher
## Created 6/5/2019
###############################################

rm(list=ls())

library("phyloseq")
library("cowplot")
library("ggthemes")
library("vegan")
library("ggplot2")

# Functions
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



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
object <- metaMDS(qd_wu, k=2, binary = FALSE)
groups_wu <- meta$nutrient

## Weighted Unifrac ehull
mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims) 
points(mds.fig, "sites", pch = 19, col = "#000000", select = control) 
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = nitrate) 
points(mds.fig, "sites", pch = 19, col = "#56B4E9", select = ammon) 
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#000000"), lwd = 2, show.groups = "Control") 
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#E69F00"), lwd = 2, show.groups = "NO3-") 
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#56B4E9"), lwd = 2, show.groups = "NH4+")
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#000000"), lwd = 2, show.groups = "Control", lty = 3) 
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#E69F00"), lwd = 2, show.groups = "NO3-", lty = 3)  
ordiellipse(object, groups_wc, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#56B4E9"), lwd = 2, show.groups = "NH4+", lty = 3) 
levels(classes) <- c("Control", "Ammonium","Nitrate")
legend("bottomrigh", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("#000000", "#56B4E9", "#E69F00"),
       y.intersp = 2, col = c("#000000", "#56B4E9", "#E69F00"))

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
levels(classes) <- c("Control","Ammonium", "Nitrate")
legend("bottomrigh", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("#000000", "#56B4E9", "#E69F00"),
       y.intersp = 2, col = c("#000000", "#56B4E9", "#E69F00"))

# Prepare data for plotting
meta <- as.data.frame(sample_data(qd))
table <- as.data.frame(otu_table(qd))
head(meta)

high <- rownames(meta[which(meta[,4] == "high"),])
control <- rownames(meta[which(meta[,4] == "Control"),])


control.df <- table[,control]
high.df <- table[,high]

test.df <- cbind(control.df, high.df)
classes <- c(rep("Control", length(control)), rep("High", length(high)))
groups <- classes
df <- test.df
dims <- c(1,2)
ellp.kind <- "ehull"

# For Weighted Unifra
qd_bj <- phyloseq::distance(qd, method = "jaccard", binary = TRUE)
object <- metaMDS(t(df), k=3, distance = "jaccard", binary = TRUE)

## Weighted Unifrac ehull
mds.fig <- ordiplot(object, display = "sites", type = "none", choices = dims) 
points(mds.fig, "sites", pch = 19, col = "#000000", select = control) 
points(mds.fig, "sites", pch = 19, col = "#E69F00", select = high) 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#000000"), lwd = 2, show.groups = "Control") 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind=ellp.kind, col = c("#E69F00"), lwd = 2, show.groups = "High") 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#000000"), lwd = 2, show.groups = "Control", lty = 3) 
ordiellipse(object, groups, conf = 0.95, label = FALSE, choices=dims, kind="se", col = c("#E69F00"), lwd = 2, show.groups = "High", lty = 3)  
levels(classes) <- c("26C","29C")
legend("bottomright", legend = levels(classes), bty = "n", pch = 21, pt.bg = c("#000000", "#E69F00"), y.intersp = 2)
