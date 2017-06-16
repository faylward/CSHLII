

#########################################################################################################################################################################
################################################################### Viral contig correlations ###########################################################################
#########################################################################################################################################################################
library(dynamicTreeCut)
library(WGCNA)
library(weights)
library(naturalsort)

# Read in relative abundance info and make one combined file
coverage_0.03 <- read.table("scaffold_full_abundance_viral_fraction_l45_p95.txt", sep="\t", header=T, quote=, row.names=1)
coverage_0.2 <- read.table("scaffold_full_abundance_cellular_fraction_l45_p95.txt", sep="\t", header=T, quote=, row.names=1)
coverage <- cbind(coverage_0.03, coverage_0.2)

names <- list()
for(i in 1:length(clust$labels)) { names[i] <- splist[i][[1]][1] }
names <- as.character(names)

# Read in TNF data
TNF <- read.table("final_phage_scaffolds.tetra", sep="\t", header=T, quote=, row.names=1)

# For weighting TNF and depth equally
TNFw <- seq(from=0.5/136, to=0.5/136, length.out=88)
depthw <- seq(from=0.5/88, to=0.5/88, length.out=136)
weights = c(depthw, TNFw)

#combine depth and tnf data into one table
dTNF <- t(cbind(coverage, TNF))

## For weighted Pearson's Correlation
corr <- wtd.cors(dTNF, weight=weights)
dist <- as.dist(1-corr)

clust <- hclust(dist, method="complete")
tree <- cutreeDynamic(clust, method="tree", cutHeight=0.005, minClusterSize=2)
colors <- labels2colors(tree)

names <- list()
for(i in 1:length(clust$labels)) { names[i] <- splist[i][[1]][1] }
names <- as.character(names)

table <- data.frame(cbind(colors, as.character(tree)))
row.names(table) <- names
colnames(table) <- c("Colors", "Clusters")
table2 <- table[order(row.names(table)),]

#Plot dendrogram
pdf(file="scaffold_cluster.pdf", height=8, width=40)
plotDendroAndColors(clust, colors, dendroLabels = names, hang = 0.05, addGuide = TRUE, guideHang = 0.05, cex.dendroLabels=0.6)
dev.off()














