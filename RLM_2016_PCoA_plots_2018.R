library("phyloseq")
library("cowplot")
library("ggthemes")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq", version = "3.8")

load(file = "qd.RData") # load the .RData phyloseq object, named 'qd'

# calculate distances from qd object
qd_un <- phyloseq::distance(qd, method = "unifrac")
qd_bj <- phyloseq::distance(qd, method = "jaccard", binary =TRUE)

# ordinate the distances
pcoa_un <- ordinate(qd, "PCoA", distance = "qd_un")
pcoa_bj <- ordinate(qd, "PCoA", distance = "qd_bj")

# plot ordination
a <- plot_ordination(qd, pcoa_bj, 
                     type = "samples", 
                     color = "nutrient", 
                     title = "Binary Jaccard")+
                     geom_point(shape = 19)+
                     scale_colour_manual(name = "Nutrient", 
                                         labels = c("Control", "Ammonium", "Nitrate"),
                                         values = c("#009292", "#490092", "#920000"))+
                     theme_few(base_size = 12)
a
b <- plot_ordination(qd, pcoa_un, 
                     type = "samples", 
                     color = "nutrient", 
                     title = "Unweighted Unifrac")+
                     geom_point(shape = 19)+
                     scale_colour_manual(name = "Nutrient", 
                                         values = c("#009292", "#490092", "#920000"))+
                     theme_few(base_size = 12)
b
plot <- plot_grid(a + theme(legend.position = "none"),
                  b + theme(legend.position = "none"), 
                  nrow = 1, align = 'w',labels = c("A","B"))
plot

legend <- get_legend(a)

final.plot <- plot_grid(plot, legend, rel_widths = c(3, .7))
final.plot

##############################################################################################################################
##### MANUSCRIPT PLOT
##############################################################################################################################

ggsave("PCOA.png", final.plot, path = "~/Documents/Dropbox/1. Research/PhD/Coral_2016/R/Figures/", width = 6, height = 3, units = "in")

