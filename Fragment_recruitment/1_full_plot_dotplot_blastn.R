
library(ggplot2)
library(scales)

#########################################################################################
###### New 25m with only G1 and all others as scaffolds- 15Kb cutoff for inclusion ######
#########################################################################################

scaff <- read.table("scaffold_index2.txt", row.names=1, sep="\t")
names <- as.character(scaff$V8)

# plot by group
#x <- read.table(file="0025m_long_output_30_besthit_norm_index_g1.txt", header=T, sep="\t")
x <- read.table(file="0025m_blastn_summary.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0025m_dotplot_blastn.pdf", height=9, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(range=c(0.5, 5)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(75, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=4, color="black"), panel.grid.major = element_blank(), panel.grid.minor.x = element_line(colour = "grey90"), panel.grid.major.y = element_blank()) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits=c(0,108), breaks=c(1:107), labels=names, expand = c(0, 0))
dev.off()







###################################
############ OLD stuff ############
###################################






################################################################
scaff <- read.table("scaffold_index.txt", row.names=1, sep="\t")
names <- row.names(scaff)
#################
###### 25m ######
#################

# plot by group
x <- read.table(file="0025m_long_output_30_besthit_norm_index.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0025m_dotplot_30_besthit_index.pdf", height=9, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_text(size=4, color="black"), panel.grid.major = element_blank(), panel.grid.minor.x = element_line(colour = "grey90"), panel.grid.major.y = element_blank()) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits=c(0,130), breaks=c(1:129), labels=names, expand = c(0, 0))
dev.off()



x <- read.table(file="0025m_long_output_30_besthit_norm.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0025m_dotplot_30_besthit_index.pdf", height=11, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90")) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits = c(0,177), expand = c(0, 0))
dev.off()

################
##### 75m ######
################
x <- read.table(file="0075m_long_output_30_besthit_norm.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))


pdf(file="0075m_dotplot_30_besthit.pdf", height=11, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90")) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits = c(0,177), expand = c(0, 0))
dev.off()

################
##### 125m #####
################
x <- read.table(file="0125m_long_output_30_besthit_norm.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0125m_dotplot_30_besthit.pdf", height=11, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90")) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits = c(0,177), expand = c(0, 0))
dev.off()

#################
##### 200m ######
#################
x <- read.table(file="0200m_long_output_30_besthit_norm.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0200m_dotplot_30_besthit.pdf", height=11, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90")) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits = c(0,177), expand = c(0, 0))
dev.off()

#################
##### 500m ######
#################
x <- read.table(file="0500m_long_output_30_besthit_norm.txt", header=T, sep="\t")
to_plot <- x[,c('cruise_index', 'scaffold_index', 'num_hits', 'percid')]
to_plot <- data.frame(as.matrix(to_plot))
print(max(to_plot$num_hits))

pdf(file="0500m_dotplot_30_besthit.pdf", height=11, width=4)
ggplot(to_plot, aes(x = cruise_index, y = scaffold_index, size=num_hits, fill=percid)) + scale_size_continuous(limits=c(0, 69000)) + geom_point(shape = 21, stroke=0, alpha=0.8) + scale_fill_gradientn(colours = c("firebrick4", "firebrick1", "grey95", "dodgerblue", "dodgerblue4"), limits=c(30, 100),  oob=squish)  + theme(axis.text.x = element_text(face="bold", color="#993333", size=8, angle=0), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.major = element_line(colour = "grey90"), panel.grid.minor = element_line(colour = "grey90")) + scale_x_continuous(breaks = c(1:12), labels=c("224", "225", "226", "227", "229", "231", "232", "233", "234", "236", "237", "238")) + scale_y_continuous(limits = c(0,177), expand = c(0, 0))
dev.off()

####################################
####################################
####################################

