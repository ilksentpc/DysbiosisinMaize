### R script process

library(phyloseq)
library(file2meco)
library(microeco)
library(MicrobiotaProcess)
library(ggplot2)
library(vegan)


#### Import data ####

setwd("~/Documents/PhD/chapter01/bacteria/leaf")

otu_leafbac <- read.csv("otu_leafbac.csv", row.names = 1)
tax_lefse <- read.csv("tax_lefseorder.csv", row.names = 1)
otumat <- as.matrix(otu_leafbac)
taxmat <- as.matrix(tax_lefse)
metadata <- read.csv("metadata_leafbac.csv", row.names = 1)
sampledata = sample_data(data.frame(metadata))
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

OTU_TREE = read_tree("bacterial_leaf.tre")
physeq = phyloseq(OTU, TAX, sampledata, OTU_TREE)
physeq = phyloseq(OTU, TAX, sampledata, OTU_TREE)
physeq

meco_dataset <- phyloseq2meco(physeq)
meco_dataset

### Differential abundance analysis-LEfSe ###

t1 <- trans_diff$new(dataset = meco_dataset, method = "lefse", group = "Genotype", alpha = 0.05, lefse_subgroup = NULL)
png("Threshold4_all.png", width=12, height=8, units="in", res=1200)
t1$plot_diff_bar(threshold = 4, group_order = c("Perennial teosinte", "B. teosinte", "MX Landrace", "MX Inbred", "US Landrace", "US Inbred"), 
                 color_values = c("Perennial teosinte"= "#1f77b4", "B. teosinte"= "#2ca02c", "MX Landrace"="#d62728", "US Landrace"="#9467bd", "MX Inbred"= "yellow", "US Inbred"="#ff7f0e"), axis_text_y = 17) +  
theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) + theme(axis.text.y = element_text(colour = "black")) +
theme(legend.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 12))
dev.off() 

### Venn Diagram

dataset1 <- meco_dataset$merge_samples(use_group = "Genotype")

t1 <- trans_venn$new(dataset1, ratio = NULL)

VennOTU <- t1$otu_table
write.csv(VennOTU, "vennOTU.csv", row.names=TRUE)

####### Community composition analysis ########

meco_dataset$sample_table$Genotype  %<>% factor(levels = c("Perennial teosinte", "B. teosinte", "MX Landrace", "MX Inbred", "US Landrace", "US Inbred"))

png("phylum_bar_leafbac.png", width=17, height=9, units="in", res=600)
t1$plot_bar(others_color = "grey70", facet = "Genotype", xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = NULL, strip_text = 7) + 
theme(axis.text.y = element_text(colour = "black", face = "bold")) + theme(legend.text = element_text(size = 15)) + 
theme(strip.text = element_text(color = "black", face = "bold"))
dev.off() 

library(dplyr)

taxonomy_totals <- t1$data_abund %>%
  group_by(Taxonomy) %>%
  summarise(Total_Abundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(Percentage = Total_Abundance / sum(Total_Abundance) * 100)
View(taxonomy_totals)

df_tax <- taxonomy_totals
write.csv(df_tax, "df_tax.csv", row.names=TRUE)

t2 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 40)

png("heatmapleafbac.png", width=17, height=9, units="in", res=600)
t2$plot_heatmap(facet = "Genotype", xtext_keep = FALSE, withmargin = FALSE, strip_text = 7, ytext_size = 14) + 
theme(axis.text.y = element_text(colour = "black")) + theme(legend.text = element_text(size = 15)) + 
theme(strip.text = element_text(color = "black"))
dev.off()


#### Weighted Unifrac Distance ###

meco_dataset$cal_betadiv(unifrac = TRUE)
meco_dataset$beta_diversity$wei_unifrac
library(magrittr)
meco_dataset$sample_table$Genotype  %<>% factor(levels = c("Perennial teosinte", "B. teosinte", "MX Landrace", "MX Inbred", "US Landrace", "US Inbred"))
t1 <- trans_beta$new(dataset = meco_dataset, group = "Genotype", measure = "wei_unifrac")
t1$cal_ordination(ordination = "PCoA")
png("wunifrac_leaf_bac2.png", width=13, height=9, units="in", res=600)
x1 <- t1$plot_ordination(plot_color = "Genotype", plot_shape = "Genotype", plot_type = c("point", "ellipse"))
x2 <- x1 + ggtitle("PCoA of Weighted UniFrac") + theme(
  text = element_text(size = 16, color = "black"),  
  axis.title = element_text(size = 18),  
  axis.text = element_text(size = 14) 
)
x2
dev.off()

#### Non-Metric Multidimensional Scaling (NMDS) ####

library(magrittr)
meco_dataset$sample_table$Genotype  %<>% factor(levels = c("Perennial teosinte", "B. teosinte", "MX Landrace", "MX Inbred", "US Landrace", "US Inbred"))
meco_dataset$cal_betadiv(method = "bray", unifrac = FALSE, binary = FALSE)
t1 <- trans_beta$new(dataset = meco_dataset, group = "Genotype", measure = "bray")

t1$cal_ordination(ordination = "NMDS", ncomp = 3)
class(t1$res_ordination)
png("nmds_leaf.png", width=13, height=9, units="in", res=600)
 z <- t1$plot_ordination(plot_color = "Genotype", plot_shape = "Genotype", plot_type = c("point", "ellipse"), 
                         NMDS_stress_pos = c(0.8, -2), NMDS_stress_text_prefix = "Stress: ")
z1 <- z + theme(
  text = element_text(size = 16, color = "black"),  
  axis.title = element_text(size = 18),  
  axis.text = element_text(size = 14) 
)
z1
dev.off()

#### Cluster dendrogram ###
library(readr)
otu_leafbac <- read.csv("otu_leafbac.csv", row.names = 1)
str(otu_leafbac)
metadata <- read.csv("metadata_leafbac.csv", row.names = 1)

library(vegan)
bray_dist <- vegdist(t(otu_leafbac), method = "bray") ##complete#muchbetter
clustOut <- hclust(bray_dist)

library(ape)
clustPhy <- as.phylo.hclust(clustOut)

png("braycurtis_cluster.png", width=9, height=8, units="in", res=1200)
p = plot.phylo(clustPhy, show.tip.label = TRUE, cex = 0.5, label.offset = 0.01)
tiplabels(tip = which(metadata$Genotype == "B. teosinte"), pch = 16, col = "orange", cex = 0.6, offset = 0.005)
tiplabels(tip = which(metadata$Genotype == "MX Inbred"), pch = 16, col = "deeppink", cex = 0.6, offset = 0.005)
tiplabels(tip = which(metadata$Genotype == "Perennial teosinte"), pch = 16, col = "forestgreen", cex = 0.5, offset = 0.005)
tiplabels(tip = which(metadata$Genotype == "US Landrace"), pch = 16, col = "green", cex = 0.6, offset = 0.005)
tiplabels(tip = which(metadata$Genotype == "MX Landrace"), pch = 16, col = "purple", cex = 0.6, offset = 0.005)
tiplabels(tip = which(metadata$Genotype == "US Inbred"), pch = 16, col = "yellow", cex = 0.6, offset = 0.005)
axisPhylo(side = 1)
mtext(side = 1, line = 2, adj = 0.5, text = "Cluster Height")
legend("topleft", legend = c("Perennial teosinte", "B. teosinte", "MX Landrace", "MX Inbred", "US Landrace", "US Inbred" ), 
       pch = 16, bty = "n", cex = 0.8, col = c("forestgreen", "orange", "purple", "deeppink", "green", "yellow", title = "Plant Genotype"))
dev.off() 






