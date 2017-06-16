library(naturalsort)
library(ggplot2)

data1 <- read.table(file="longform_data_newdata_for_comparison.txt", sep="\t", header=T)
data1$metagenome <- log10(data1$metagenome * 100)
data1$virome <- log10(data1$virome * 100)

pdf(file="Scaffold_mg_vs_vi_bubbleplot.pdf", height=8, width=8)
ggplot(data1, aes(x = virome, y = metagenome, size = size, alpha=0.7)) + geom_point(shape = 21, colour = "blue", fill = "blue") + ggtitle("Viral Abundance: Metagenome vs Virome") + labs(x = "log10 Sum Virome Coverage", y = "Log10 Sum Metagenome Coverage") + geom_abline(intercept = 0) + xlim(1, 7.5) + ylim(1, 6.5) + theme(axis.text.x = element_text(color="black", size=12, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=12, color="black"), panel.grid.major = element_blank(), panel.grid.major.x = element_line(colour = "grey70"), panel.grid.major.y = element_line(colour = "grey70"), panel.border = element_rect(colour = "grey70", fill=NA, size=1))
dev.off()
#################
#################
#################
