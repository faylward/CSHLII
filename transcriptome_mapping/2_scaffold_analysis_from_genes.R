
library(naturalsort)
library(gplots)
library(ggplot2)
########################################################################################################
################ Load data, get scaffold readcounts from gene counts, and normalize ####################
########################################################################################################
#raw <- read.table(file="genes_80pctid_last_idx.txt", header=T, row.names=1, sep="\t")
raw <- read.table(file="genes_90pctid_last_idx.txt", header=T, row.names=1, sep="\t")
#raw <- read.table(file="genes_95pctid_last_idx.txt", header=T, row.names=1, sep="\t")
raw <- raw[naturalsort(row.names(raw)),]

to_remove <- c("scaffold399|size4924_1")  #, "scaffold50|size23177_9", "scaffold411|size8398_1")
to_include <- setdiff(row.names(raw), to_remove)
counts <- raw[to_include,]

counts1 <- counts[,4:47]
ag <- aggregate(counts1, by=list(counts$scaffold_idx), FUN=sum)

lkp <- read.table(file="scaffold2group.txt", header=T, row.names=1)
lookup <- lkp[naturalsort(row.names(lkp)),]
#scaff <- intersect(row.names(lookup), row.names(counts))
#lookup2 <- lookup[scaff,]

coords <- match(ag$Group.1, lookup$group_idx)
coords2 <- coords[!is.na(coords)]
names <- as.character(lookup$Final_Group)[coords2]
row.names(ag) <- names
fcounts <- ag[,2:45]

# sort scaffolds and filter out low-abundance ones
rowsum <- rowSums(fcounts)
pdtable <- cbind(fcounts, rowsum)
pdtable1 <- pdtable[order(pdtable$rowsum, decreasing=T),]
abund <- row.names(pdtable1)[which(pdtable1$rowsum >= 44)]
abund2 <- row.names(pdtable1)[which(pdtable1$rowsum >= 88)]
excluded <- setdiff(lookup$Final_Group, abund)
pdtable2 <- pdtable1[abund, ]
ntable <- pdtable2[,1:44]

# normalize transcripts with molecular standard measurements
stan <- read.table(file='standard_normalization.info2', header=TRUE, row.names=1, sep="\t", quote="")
norm_factor <- as.numeric(seq(from=1, to=1, by=1)) / stan$SNC
norm_table <- data.frame(scale(ntable, scale=norm_factor, center=F))



########################################################################
################# Plot total transcripts/L of viruses ##################
########################################################################
stan2 <- read.table(file='sample_times.txt', header=TRUE, row.names=1, sep="\t", quote="")
times <- stan2$Time3
lcounts <- fcounts * 100000 
colsum <- colSums(lcounts)

to_plot <- data.frame(cbind(times, colsum))
colnames(to_plot) <- c("times", "transcripts")
to_plot$transcripts <- log10(to_plot$transcripts)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}
rdata <- summarySE(to_plot, measurevar="transcripts", groupvars="times")

pdf(file="Viral_transcripts_per_liter_line_90pctid.pdf", height=3.5, width=7)
a <- geom_errorbar(data=rdata, aes(x=times, ymin=transcripts-se, ymax=transcripts+se), width=0.7, size=1, color="firebrick")
b <- geom_line(data = rdata, aes(x=times, y=transcripts), color="firebrick") 
c <- geom_point(data=rdata, aes(x=times, y=transcripts), color="firebrick", size=3)
d <- theme(axis.text.x=element_text(color="#993333", hjust=1, size=12), axis.text.y=element_text(color="black", hjust=1, size=7), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color="grey90"), panel.border = element_rect(colour = "grey80", fill=NA, size=1))
e <- scale_x_continuous(limits=c(0, 24), breaks=rdata$times, labels=c('0200', '0600', '1000', '1400', '1800', '2200'))
ggplot() + a + b + c + d + e
dev.off()

subset_a <- subset(to_plot, to_plot$times %in% c(2, 6, 10))
subset_b <- subset(to_plot, to_plot$times %in% c(14, 18, 22))

wr <- wilcox.test(subset_a$transcripts, subset_b$transcripts, alternative="l")

#t.test(subset_a$transcripts, subset_b$transcripts)


##############################################################
################# get full normalized table ##################
##############################################################
ntable2 <- pdtable1[,1:44]
full_norm_table <- data.frame(scale(ntable2, scale=norm_factor, center=F))
means <- apply(ntable2, 1, mean)
full_norm_table2 <- cbind(means, full_norm_table)

excluded2 <- setdiff(lookup$Final_Group, row.names(full_norm_table2))
empty <- matrix(0, nrow=length(excluded2), ncol=dim(full_norm_table2)[2])
colnames(empty) <- colnames(full_norm_table2); row.names(empty) <- excluded2
all <- rbind(full_norm_table2, empty)
write.table(all, "viral_scaffolds_normcounts.tbl", sep="\t", quote=F)
write.table(all[,1], "viral_scaffolds_mean_normcounts.tbl", sep="\t", quote=F, row.names=row.names(all))

################################################################################
############### Run RAIN Analysis to test for 24-hr periodicity ################
################################################################################
library(rain)
######### Get a regularized T/F time series vector for RAIN ###############
t <- stan$Time
ft <- seq(from=4, to=208, by=4)

for (i in 1:length(ft)) { 
	if (ft[i] %in% t) {
		ft[i] = 1
	}
	else {
		ft[i] = 0
	}
}
 
# Run on all scaffolds
stan_final <- rain(t(norm_table), period=24, measure.sequence=ft, deltat=4)
fdr_final <- p.adjust(as.numeric(stan_final$pVal), method="fdr")	
signif_final <- sum(fdr_final <= 0.1)
print(signif_final)

###############################################
############ Get peak time ####################
###############################################

## Define Function
gene.glm.periodic <- function (gene.counts) {

	Result <- list()
	sample.times <- as.numeric(c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204, 208))
	gene.test <- lm(as.numeric(gene.counts) ~ cos(2*pi*sample.times/24) + sin(2*pi*sample.times/24))
	LR.obs <- gene.test$null.deviance - gene.test$deviance
	
	coef <- gene.test$coefficients
	amp <- sqrt( coef[2]^2 + coef[3]^2 )
	phase <- atan( -1 * coef[3] / coef[2] )

	phase <- atan( -1 * coef[3] / coef[2] )
	if (phase > 0) { peak <- 24 - phase*24/(2*pi) }
	else {  peak <- -1*phase*24/(2*pi) }
	if (coef[2] < 0) { 
		if (peak < 12 ) {
			peak <- peak + 12
		} else {
			peak <- peak - 12
		}
	}
	Result$fit <- gene.test
	Result$amp <- amp
	Result$peak <- peak
	Result$pval <-  pchisq(LR.obs, 2, lower.tail=FALSE)
	return(Result)
}
# End

# Run function on data
data1 <- apply(norm_table, 1, gene.glm.periodic)
length(data1)
peaks <- list()
for(i in 1:length(data1)) { peaks[i] <- data1[[i]]$peak }
peaks <- as.numeric(peaks)
head(peaks)

######################################################################################
############################# print out final data tables ############################
######################################################################################
diel <- fdr_final <= 0.1
diel_table <- cbind(diel, peaks, fdr_final, stan_final)
diel_table <- diel_table[order(diel_table$fdr_final),]

other <- matrix(NA, nrow=length(excluded), ncol=dim(diel_table)[2])
colnames(other) <- colnames(diel_table); row.names(other) <- excluded
full_diel_table <- rbind(diel_table, other)
write.table(full_diel_table, "viral_scaffolds_from_genes_full_diel_table.tbl", sep="\t", quote=F)

############################################################################################
##################### Make a violin plot for most abundant scaffolds/groups ################
############################################################################################
#long_scaffolds <- as.character(unique(lookup$Final_Group[which(lookup$length > 15000)]))
ls <- as.character(unique(lookup$Final_Group[which(lookup$length > 10000)]))
long_scaffolds <- intersect(ls, abund2)
scaff <- intersect(row.names(norm_table), long_scaffolds)

# find out which of these scaffolds are diel
diel_table <- read.table("viral_scaffolds_from_genes_full_diel_table.tbl", sep="\t", header=T, row.names=1)
d <- subset(diel_table, diel_table$diel==TRUE & row.names(diel_table) %in% scaff)

thresh <- range(stan$SNC * 100000)

med <- apply(norm_table, 1, mean)
t <- data.frame(cbind(norm_table, med))
t <- t[scaff,]
t <- t * 100000
ord <- t[order(t$med, decreasing=T),]
subset <- ord#[1:83,]
subset_ord <- subset[order(subset$med), 1:44]
mean_thresh <- mean(thresh[1], thresh[2])
subset_ord[subset_ord == 0] <- thresh[1]

library(reshape2)
library(ggplot2)
d <- melt(t(subset_ord))

# plot by group
to_plot <- data.frame(cbind(as.numeric(d$Var2), d[,3]))
colnames(to_plot) <- c("index", "abund")
to_plot$abund <- log10(to_plot$abund)

n <- strsplit(row.names(subset_ord), "|", fixed=T)
names <- list()
for(i in 1:length(n)) { names[i] <- n[i][[1]][[1]] }
names <- as.character(names)
names <- gsub("scaffold", "VS", names)

pdf(file="viral_scaffolds_transcriptome_violin.pdf", height=7.5, width=4.5)
#x <- geom_rect(data=to_plot, aes(ymin=log10(thresh[1]), ymax=log10(thresh[2]), xmin=0, xmax=83), fill='grey90', colour=NA)
x <- geom_hline(aes(yintercept=log10(thresh[2])), linetype='dashed', size=0.5, colour='firebrick')
y <- geom_hline(aes(yintercept=log10(thresh[1])), linetype='dashed', size=0.5, colour='firebrick')
a <- geom_violin(data=to_plot, aes(x = index, y = abund, group = cut_width(index, 1)), fill= "dodgerblue4", color="dodgerblue4", trim=T)
b <- scale_y_continuous(expand=c(0.1, 0.1)) # limits=c(3, 8)
c <- theme(axis.text.x=element_text(color="#993333", hjust=1, size=12), axis.text.y=element_text(color="black", hjust=1, size=7), panel.background=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = "grey90"), panel.grid.major.y = element_line(colour = "grey90"), panel.border = element_rect(colour = "grey80", fill=NA, size=1))
d <- scale_x_continuous(breaks = c(1:60), labels=names, expand=c(0.01, 0.01))
e <- coord_flip()
ggplot() + x + y + a + b + c + d + e
#ggplot() + a + c + e
dev.off()

#######################################################################################
################ Plot histogram of when diel genes are peaking ########################
#######################################################################################
long_ndiel <-  diel_table$peaks[which(row.names(diel_table) %in% long_scaffolds & diel_table$fdr_final > 0.1)]
long_diel <-  diel_table$peaks[which(row.names(diel_table) %in% long_scaffolds & diel_table$fdr_final <= 0.1)]

diel_peaks <- diel_table$peaks[which(diel_table$diel==TRUE)]
ndiel_peaks <- diel_table$peaks[which(diel_table$diel==FALSE)]
pdf("diel_genes_peak_times.pdf", height=4, width=6)
hist(diel_peaks, col="dodgerblue2", breaks=24, axes=F, xlim=c(0, 24), xlab="", ylab="", main="")
#par(new=T)
#hist(long_diel, col="dodgerblue", breaks=24, axes=F)
#plot(density(diel_peaks), col="dodgerblue", type="l")
#abline(v=9.25)
axis(1, at=c(0, 4, 8, 12, 16, 20, 24), labels=c("0000", "0400", "0800", "1200", "1600", "2000", "2400"))
axis(2, las=2)
box(bty="l")
dev.off()
#######################
#######################
#######################



#######################################################################################
############################ Plot diel dotplot ########################################
#######################################################################################

diel <- cbind(diel_peaks, rep(1, length(diel_peaks)), rep("diel", length(diel_peaks)))
ndiel <- cbind(ndiel_peaks, rep(1, length(ndiel_peaks)), rep("irreg", length(ndiel_peaks)))
full_data <- data.frame(rbind(diel, ndiel))
colnames(full_data) <- c("peak_time", "index", "diel")
full_data$peak_time <- as.numeric(as.character(full_data$peak_time))
full_data$index <- as.numeric(as.character(full_data$index))

a <- geom_rect(data=full_data, aes(xmin=-0.5, xmax=6, ymin=0.7, ymax=1.3), fill='grey92', colour=NA)
b <- geom_rect(data=full_data, aes(xmin=18.5, xmax=24.5, ymin=0.7, ymax=1.3), fill='grey92', colour=NA)

c <- geom_point(data=full_data[which(full_data$diel == "irreg"),], position = position_jitter(w = 0, h = 0.15), aes(peak_time, index), alpha=0.45, shape=16, size=3, colour="grey60")
d <- geom_point(data=full_data[which(full_data$diel == "diel"),], position = position_jitter(w = 0, h = 0.15), aes(peak_time, index), alpha=0.85, shape=16, size=4, colour="dodgerblue2")

e <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey20', size=16), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
f <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
g <- scale_y_continuous(breaks=c(1), labels="all_scaffolds", expand = c(0.01, 0))
h <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
i <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

pdf("diel_scaffolds_dotplot.pdf", height=1, width=8)
ggplot() + a + b + c + d + e + f + g + h + i
dev.off()










