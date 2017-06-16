

library(naturalsort)
library(reshape2)
library(ggplot2)
library(scales)
library(data.table)
#########################################################################################################################################################################################
######################################################### Viral scaffold normalized average coverage profiles ###########################################################################
#########################################################################################################################################################################################

# Read in relative abundance info and make one combined file
coverage_0.03 <- read.table("scaffold_rawcov_0.03_fltr_l45_p95.txt", sep="\t", header=T, quote=, row.names=1)
coverage_0.2  <- read.table("scaffold_rawcov_0.2_fltr_l45_p95.txt", sep="\t", header=T, quote=, row.names=1)

total_vmbp <- data.frame(read.table("virome_readpool.bp", sep="\t", header=T, row.names=1))
tv <- total_vmbp[naturalsort(row.names(total_vmbp)), , drop=F]

total_cmbp <- read.table("metagenome_readpool.bp", sep="\t", header=T, row.names=1)
tc <- total_cmbp[naturalsort(row.names(total_cmbp)), , drop=F]

cellcov <- scale(coverage_0.2, scale=tc$bp/10000000000, center=F)
virocov <- scale(coverage_0.03, scale=tv$bp/10000000000, center=F)

ca <- cellcov + virocov
ca2 <- ca[naturalsort(row.names(ca)),]

ratio <- virocov/cellcov
ratio2 <-  ratio[naturalsort(row.names(ratio)),]
ratio2[ratio2==NaN] <- NA
ratio2[ratio2==Inf] <- NA

### Print full data tables for SI material
write.table(ratio2, "scaffold_vc_ratio.txt", sep="\t", quote=F)
write.table(virocov, file="scaffold_full_abundance_viral_fraction.txt", sep="\t", quote=F)
write.table(cellcov, file="scaffold_full_abundance_cellular_fraction.txt", sep="\t", quote=F)

# print mean abundance and ratio values for SI material
vir <- apply(virocov, 1, mean)
cel <- apply(cellcov, 1, mean)
r <- apply(ratio2, 1, FUN=mean, na.rm=T)
write.table(vir, file="scaffold_mean_abundance_viral_fraction.txt", sep="\t", quote=F)
write.table(cel, file="scaffold_mean_abundance_cellular_fraction.txt", sep="\t", quote=F)
write.table(r, file="scaffold_mean_abundance_vc_ratio.txt", sep="\t", quote=F)

## insert gap in the time-series and output tables for plotting
gap <- data.frame(matrix(data = 0, nrow = dim(ca2)[1], ncol = 2))
colnames(gap) <- c("gap1", "gap2")
ca3 <- cbind(ca2[,1:25], gap, ca2[,26:44])
ca3[,"S69C001"] <- 0

maxval <- range(ratio2, finite=T)[2]
ratio2[ratio2==Inf] <- maxval

ratio3 <- cbind(ratio2[,1:25], gap, ratio2[,26:44])
ratio3[,"S69C001"] <- 0

write.table(ca3, file="scaffold_combined_abundance_for_plotting.txt", sep="\t", quote=F)
write.table(ratio3, file="scaffold_ratio_for_plotting.txt", sep="\t", quote=F)













