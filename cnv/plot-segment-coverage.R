options(warn=1)
library(DNAcopy)
library(optparse)

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--tumor", type="character", help="tumor coverage data file"),
		make_option("--normal", type="character", help="remission coverage data file"),
		make_option("--output", type="character", help="output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$patient)) stop("patient ID not specified")
if (is.na(opt$tumor)) stop("tumor coverage file not specified")
if (is.na(opt$normal)) stop("normal coverage file not specified")
if (is.na(opt$output)) stop("output file not specified")

# for test purposes
#opt <- data.frame(patient="1020540_dia", tumor="../reseq/cnv/segmented_coverage/1020540_Diagnosis.segmented-coverage.tsv", normal="../reseq/cnv/segmented_coverage/1020540_Remission.segmented-coverage.tsv", stringsAsFactors=F)

# normalization factor (correct for different read count)
normal.chrs <- list(
		"1020540_dia" = c("chr7", "chr8", "chr9", "chr10"),
		"NH17331_dia" = c("chr1", "chr2", "chr3"), 
		"944_dia" = c("chr1", "chr2", "chr3"), 
		"944_rel" = c("chr1", "chr2", "chr3"), 
		"HS6_dia" = c("chr2", "chr3", "chr5"), 
		"933_dia" = c("chr1", "chr2", "chr3"), 
		"933_rel" = c("chr4", "chr9", "chr16"), 
		"1019357_rel" = c("chr1", "chr2", "chr3"), 
		"1009302_dia" = c("chr1", "chr2", "chr3"), 
		"1009302_dia" = c("chr1", "chr2", "chr3"),
		"1004564_dia" = c("chr1", "chr2", "chr3"),
		"1010781_dia" = c("chr2", "chr3", "chr5"),
		"1017005_rel" = c("chr2", "chr3", "chr16"),
		"1019357_dia" = c("chr1", "chr2", "chr3"),
		"1019964_dia" = c("chr3", "chr11", "chr19"),
		"1020076_dia" = c("chr2", "chr3", "chr4"),
		"1020540_rel" = c("chr1", "chr2", "chr3"),
		"1020583_dia" = c("chr2", "chr3", "chr16"),
		"1020583_rel" = c("chr1", "chr2", "chr3"),
		"1021087_dia" = c("chr1", "chr2", "chr3"),
		"1021392_dia" = c("chr1", "chr2", "chr3"),
		"1021392_rel" = c("chr1", "chr2", "chr3"),
		"1021631_dia" = c("chr2", "chr3", "chr16"),
		"1021865_rel" = c("chr2", "chr3", "chr16"),
		"1022914_dia" = c("chr13", "chr15", "chr16"),
		"1022914_rel" = c("chr13", "chr15", "chr16"),
		"1023056_dia" = c("chr2", "chr3", "chr16"),
		"1023056_rel" = c("chr2", "chr3"),

		"331_dia" = c("chr2", "chr3", "chr16"),
		"777_dia" = c("chr1", "chr2", "chr3")) 

t <- read.delim(opt$tumor, header=F)
n <- read.delim(opt$normal, header=F)
m <- merge(n, t, by=c("V1", "V2"))
m$ratio <- log((m[,"V3.y"] / sum(m[,"V3.y"])) / (m[,"V3.x"] / sum(m[,"V3.x"])), 2)
m$ratio <- m$ratio - ifelse(!is.null(normal.chrs[[opt$patient]]), mean(m$ratio[which(m$V1 %in% normal.chrs[[opt$patient]] & is.finite(m$ratio))]), 0)

set.seed(25)
CNA.object <-CNA(genomdat = m[,"ratio"], chrom = m[,"V1"], maploc = m[,"V2"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2)

pdf(opt$output, width=30, paper="A4r")
par(mfrow=c(5,5), mar=c(0.5,2,1.5,1))
for (c in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
{
	plot(m[m$V1==c,"V2"], m[m$V1==c,"ratio"], ylim=c(-1.5,1.5), main=c, cex=0.3, xaxt='n')
	abline(0, 0)
	for(i in which(segs$output$chrom==c)) {
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
	}
}
dev.off()
