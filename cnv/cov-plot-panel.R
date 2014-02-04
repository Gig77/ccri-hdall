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
		   "446.dia" = 0, "446.rel" = 0.8, 
		   "460.dia" = 0, "460.rel" = 0, 
		   "X.dia" = 0.4, "X.rel" = 1, 
		   "786.dia" = 0.6, "786.rel" = 1.1, 
		   "MB2.dia" = 0.8, "MB2.rel" = 0.1, 
		   "KE12025.dia" = 1.45, "KE12025.rel" = 0.35, 
		   "C.dia" = 0.05, "C.rel" = 0.2, 
		   "243.dia" = 1, "243.rel" = 1.7, 
		   "545.dia" = 0, "545.rel" = 0, 
		   "818.dia" = 0.5, "818.rel" = 0.05, 
		   "1024543.dia" = 0.25, "1024543.rel" = 0.75,
		   "592.dia" = -0.15, "592.rel" = 0.5, 
		   "933.dia" = -1.5, "933.rel" = -0.5, 
		   "944.dia" = 0.7, "944.rel" = 1.7, 
		   "1010661.dia" = 0.5, 
		   "HL1.rel" = -2.2, 
			"FE1.rel" = 1.1, 
			"GD1.rel" = 1.0, 
			"BL16.rel" = 1.0, 
			"KL16.rel" = 0.8) 
   
# contamination factor (correct for normal contamination)
cf <- list("430.dia" = 1.0, "430.rel" = 1.0, 
		   "314.dia" = 1.0, "314.rel" = 1.0, 
		   "399.dia" = 1.0, "399.rel" = 1.0, 
		   "446.dia" = 1.0, "446.rel" = 1.1, 
		   "460.dia" = 1.0, "460.rel" = 1.0, 
		   "545.dia" = 1.0, "545.rel" = 1.0, 
		   "818.dia" = 1.0, "818.rel" = 1.0, 
		   "X.dia" = 1.0, "X.rel" = 0.6, 
		   "786.dia" = 1.0, "786.rel" = 0.7, 
		   "1024543.dia" = 1.7, "1024543.rel" = 1.3,
		   "KE12025.dia" = 1.1, "KE12025.rel" = 1.2, 
		   "592.dia" = 1.0, "592.rel" = 1.0,
		   "1010661.dia" = 0.7, 
		   "LM18.rel" = 0.4, 
   		   "BM18.rel" = 0.4) 
   chr.order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

if (is.na(opt$patient)) stop("patient not specified")
if (is.null(opt$diagnosis) & is.null(opt$relapse)) stop("diagnosis or relapse sample not specified")
if (is.na(opt$remission)) stop("remission sample not specified")

rem <- read.delim(opt$remission, header=F)
names(rem) <- c("chr", "start", "end", "name", "length", "strand", "total.rem", "avg.rem")
rectangles <- data.frame(xmin=c(9, 18, 33, 43, 53, 60, 68, 72, 81, 86, 89, 96), xmax=c(14, 24, 37, 47, 56, 66, 70, 74, 83, 88, 93, 97), ymin=-1.5, ymax=1.5)

if (!is.null(opt$diagnosis)) {
	dia <- read.delim(opt$diagnosis, header=F)
	names(dia) <- c("chr", "start", "end", "name", "length", "strand", "total.dia", "avg.dia")
	m.dia <- merge(rem, dia, by=c("chr", "start", "end", "name", "length", "strand"))
	m.dia <- cbind(m.dia, sapply(strsplit(as.character(m.dia$name), "-"), "[", 1))
	#m.dia <- cbind(m.dia, (log(m.dia$avg.dia/m.dia$avg.rem, 2) + nf[[paste(opt$patient, ".dia", sep="")]]) * cf[[paste(opt$patient, ".dia", sep="")]])
	m.dia <- cbind(m.dia, log(m.dia$avg.dia/m.dia$avg.rem, 2))
	m.dia <- m.dia[with(m.dia, order(chr, start)),]
	m.dia[,12] <- m.dia[,12] + ifelse(!is.null(nf[[paste(opt$patient, ".dia", sep="")]]), nf[[paste(opt$patient, ".dia", sep="")]], 0)
	m.dia[,12] <- m.dia[,12] * ifelse(!is.null(cf[[paste(opt$patient, ".dia", sep="")]]), cf[[paste(opt$patient, ".dia", sep="")]], 1)
	names(m.dia)[11] <- "gene"
	names(m.dia)[12] <- "log.rem.dia"
	
	# reorder gene factor by chromosome coordinate
	genepos <- merge(aggregate(start ~ gene, m.dia, min), m.dia)[,c(3,2,1)]
	m.dia$gene <- factor(m.dia$gene, levels(m.dia$gene)[with(genepos, order(sapply(chr, function(x) { which(chr.order==x) }), start))])

	dia.plot <- ggplot() + scale_x_discrete() +
			geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) + 
			geom_point(data=m.dia, aes(gene, log.rem.dia), na.rm=T, alpha=0.6) + ylim(-1.5, 1.5) + 
			theme_bw() + 
			geom_hline(yintercept = 0, colour="black") + 
			geom_hline(yintercept = 0.58, colour="orange") + 
			geom_hline(yintercept = 1, colour="blue") + 
			geom_hline(yintercept = -1, colour="red") + 
			theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste(opt$patient, "dia"))
}

if (!is.null(opt$relapse)) {
	rel <- read.delim(opt$relapse, header=F)
	names(rel) <- c("chr", "start", "end", "name", "length", "strand", "total.rel", "avg.rel")
	m.rel <- merge(rem, rel, by=c("chr", "start", "end", "name", "length", "strand"))
	m.rel <- cbind(m.rel, sapply(strsplit(as.character(m.rel$name), "-"), "[", 1))
	#m.rel <- cbind(m.rel, (log(m.rel$avg.rel/m.rel$avg.rem, 2) + nf[[paste(opt$patient, ".rel", sep="")]]) * cf[[paste(opt$patient, ".rel", sep="")]])
	m.rel <- cbind(m.rel, log(m.rel$avg.rel/m.rel$avg.rem, 2))
	m.rel <- m.rel[with(m.rel, order(chr, start)),]
	m.rel[,12] <- m.rel[,12] + ifelse(!is.null(nf[[paste(opt$patient, ".rel", sep="")]]), nf[[paste(opt$patient, ".rel", sep="")]], 0)
	m.rel[,12] <- m.rel[,12] * ifelse(!is.null(cf[[paste(opt$patient, ".rel", sep="")]]), cf[[paste(opt$patient, ".rel", sep="")]], 1)
	names(m.rel)[11] <- "gene"
	names(m.rel)[12] <- "log.rem.rel"
	
	# reorder gene factor by chromosome coordinate
	genepos <- merge(aggregate(start ~ gene, m.rel, min), m.rel)[,c(3,2,1)]
	m.rel$gene <- factor(m.rel$gene, levels(m.rel$gene)[with(genepos, order(sapply(chr, function(x) { which(chr.order==x) }), start))])
		
	rel.plot <- ggplot() + scale_x_discrete() + 
			geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', alpha=0.4) + 
			geom_point(data=m.rel, aes(gene, log.rem.rel), na.rm=T, alpha=0.6) + ylim(-1.5, 1.5) + 
			theme_bw() + 
			geom_hline(yintercept = 0, colour="black") + 
			geom_hline(yintercept = 0.58, colour="orange") + 
			geom_hline(yintercept = 1, colour="blue") + 
			geom_hline(yintercept = -1, colour="red") + 
			theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste(opt$patient, "rel"))
}


pdf(paste("cov-plot.", opt$patient, ".reseq.pdf.part", sep=""), width=15, height=10)
if (!is.null(opt$diagnosis) & !is.null(opt$relapse)) {
	par(mfrow=c(2,1))
	grid.arrange(dia.plot, rel.plot)
} else if (!is.null(opt$diagnosis)) {
	dia.plot
} else {
	rel.plot
}
dev.off()

