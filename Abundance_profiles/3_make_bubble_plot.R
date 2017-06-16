
library(ggplot2)
library(scales)
data <- read.table(file="longform_data_newdata.txt", header=T, row.names=1) 
data$abundance[which(data$abundance > 1000)] <- 1000

names <- read.table(file="plotting_index2.txt", header=T, row.names=1, sep="\t")

# use log10 VC ratio
data$vc <- log10(data$vc)
data$vc[data$vc == -Inf] <- NA
########################################################################################################
######################### make bubbleplot using a 1200 abundance top ###################################
########################################################################################################
# with x-axis breaks at noon
pdf(file="scaffold_abundance_log10.pdf", height=9, width=4)
ggplot(data, aes(x = cast_index, y = group_index, size=abundance, fill=vc)) + scale_size_continuous(range=c(0, 4), limits=c(0, 1200)) + geom_point(shape = 21, stroke=0.2, alpha=0.7) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey90", "lightskyblue2", "lightskyblue", "dodgerblue", "dodgerblue2", "dodgerblue3", "dodgerblue4"))  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=4, color="black"), panel.grid.major = element_blank(), panel.grid.major.x = element_line(colour = "grey90"), panel.grid.major.y = element_blank()) + scale_y_continuous(expand = c(0.01, 0), labels=rev(row.names(names)), breaks=c(1:length(row.names(names)))) + scale_x_continuous(breaks=c(2.5, 8.5, 14.5, 20.5, 32.5, 38.5, 44.5))
dev.off()
#########
