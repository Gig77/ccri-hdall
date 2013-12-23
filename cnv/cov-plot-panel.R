options(warn=1)
library(ggplot2)
library(optparse)
library(gridExtra)

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--diagnosis", type="character", help="diagnosis coverage data file"),
		make_option("--relapse", type="character", help="relapse coverage data file"),
		make_option("--remission", type="character", help="remission coverage data file")
)
opt <- parse_args(OptionParser(option_list=option_list))

# normalization factor (correct for different read count)
nf <- list("430.dia" = 0.2, "430.rel" = 0.5, 
		   "314.dia" = 0, "314.rel" = 0, 
		   "399.dia" = 0, "399.rel" = 0, 
		   "446.dia" = 0, "446.rel" = 0, 
		   "460.dia" = 0, "460.rel" = 0, 
		   "545.dia" = 0, "545.rel" = 0, 
		   "818.dia" = 0.5, "818.rel" = 0.05, 
		   "1024543.dia" = 0.25, "1024543.rel" = 0.75,
		   "592.dia" = -0.15, "592.rel" = 0.5) 

# contamination factor (correct for normal contamination)
cf <- list("430.dia" = 1.0, "430.rel" = 1.0, 
		   "314.dia" = 1.0, "314.rel" = 1.0, 
		   "399.dia" = 1.0, "399.rel" = 1.0, 
		   "446.dia" = 1.0, "446.rel" = 1.0, 
		   "460.dia" = 1.0, "460.rel" = 1.0, 
		   "545.dia" = 1.0, "545.rel" = 1.0, 
		   "818.dia" = 1.0, "818.rel" = 1.0, 
		   "1024543.dia" = 1.7, "1024543.rel" = 1.3,
   		   "592.dia" = 1.0, "592.rel" = 1.0) 
chr.order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

if (is.na(opt$patient)) stop("patient not specified")
if (is.null(opt$diagnosis)) stop("diagnosis sample not specified")
if (is.na(opt$relapse)) stop("relapse sample not specified")
if (is.na(opt$remission)) stop("remission sample not specified")

dia <- read.delim(opt$diagnosis, header=F)
rel <- read.delim(opt$relapse, header=F)
rem <- read.delim(opt$remission, header=F)

names(dia) <- c("chr", "start", "end", "name", "length", "strand", "total.dia", "avg.dia")
names(rem) <- c("chr", "start", "end", "name", "length", "strand", "total.rem", "avg.rem")
names(rel) <- c("chr", "start", "end", "name", "length", "strand", "total.rel", "avg.rel")

m <- merge(rem, dia, by=c("chr", "start", "end", "name", "length", "strand"))
m <- merge(m, rel, by=c("chr", "start", "end", "name", "length", "strand"))
m <- cbind(m, sapply(strsplit(as.character(m$name), "-"), "[", 1))
m <- cbind(m, (log(m$avg.dia/m$avg.rem, 2) + nf[[paste(opt$patient, ".dia", sep="")]]) * cf[[paste(opt$patient, ".dia", sep="")]])
m <- cbind(m, (log(m$avg.rel/m$avg.rem, 2) + nf[[paste(opt$patient, ".rel", sep="")]]) * cf[[paste(opt$patient, ".rel", sep="")]])
m <- m[with(m, order(chr, start)),]

names(m)[13] <- "gene"
names(m)[14] <- "log.rem.dia"
names(m)[15] <- "log.rem.rel"

# reorder gene factor by chromosome coordinate
genepos <- merge(aggregate(start ~ gene, m, min), m)[,c(3,2,1)]
m$gene <- factor(m$gene, levels(m$gene)[with(genepos, order(sapply(chr, function(x) { which(chr.order==x) }), start))])

rectangles <- data.frame(xmin=c(9, 18, 33, 43, 53, 60, 68, 72, 81, 86, 89, 96), xmax=c(14, 24, 37, 47, 56, 66, 70, 74, 83, 88, 93, 97), ymin=-1.5, ymax=1.5)
#rectangles <- data.frame(xmin=m$gene, xmax=m$gene, ymin=-1.5, ymax=1.5)

pdf(paste("cov-plot.", opt$patient, ".reseq.pdf.part", sep=""), width=15, height=10)
par(mfrow=c(2,1))
#ggplot() + geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.8) + geom_point(m, aes(gene, log.rem.dia)) + ylim(-1.5, 1.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dia.plot <- ggplot() + scale_x_discrete() + 
		geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) + 
		geom_point(data=m, aes(gene, log.rem.dia), na.rm=T, alpha=0.6) + ylim(-1.5, 1.5) + 
		theme_bw() + 
		geom_hline(yintercept = 0, colour="black") + 
		geom_hline(yintercept = 0.58, colour="orange") + 
		geom_hline(yintercept = 1, colour="blue") + 
		geom_hline(yintercept = -1, colour="red") + 
		theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + 
		ggtitle(paste(opt$patient, "dia"))

rel.plot <- ggplot() + scale_x_discrete() + 
		geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) + 
		geom_point(data=m, aes(gene, log.rem.rel), na.rm=T, alpha=0.6) + ylim(-1.5, 1.5) + 
		theme_bw() + 
		geom_hline(yintercept = 0, colour="black") + 
		geom_hline(yintercept = 0.58, colour="orange") + 
		geom_hline(yintercept = 1, colour="blue") + 
		geom_hline(yintercept = -1, colour="red") + 
		theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + 
		ggtitle(paste(opt$patient, "rel"))
grid.arrange(dia.plot, rel.plot)
dev.off()

