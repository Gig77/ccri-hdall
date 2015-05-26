options(warn=1)
library(DNAcopy)
library(optparse)

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--case", type="character", help="coverage data file case patient (remission sample)"),
		make_option("--control", type="character", help="coverage data file control patient (remission sample, most NOT be down patient"),
		make_option("--output", type="character", help="output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- list(patient="839C", case="segmented_coverage/839C.segmented-coverage.tsv", control="segmented_coverage/92C.segmented-coverage.tsv", output="segmented_coverage/test.pdf")

if (is.na(opt$patient)) stop("patient ID not specified")
if (is.na(opt$case)) stop("case coverage file not specified")
if (is.na(opt$control)) stop("control coverage file not specified")
if (is.na(opt$output)) stop("output file not specified")

# for test purposes
#opt <- data.frame(patient="1020540_dia", tumor="../reseq/cnv/segmented_coverage/1020540_Diagnosis.segmented-coverage.tsv", normal="../reseq/cnv/segmented_coverage/1020540_Remission.segmented-coverage.tsv", circos="/mnt/projects/hdall/results/cnv/circos/1020540_dia.cnv.circos.tsv", stringsAsFactors=F)

t <- read.delim(opt$case, header=F, colClasses=c("factor", "integer", "numeric"))
n <- read.delim(opt$control, header=F, colClasses=c("factor", "integer", "numeric"))
m <- merge(n, t, by=c("V1", "V2"))
m.chr21 <- m[m$V1 %in% c("chr21", "21"),]
m.other <- m[!(m$V1 %in% c("chr21", "21")),]
m.chr21$ratio <- log2(m.chr21[,"V3.y"] / m.chr21[,"V3.x"]) - log2(sum(m.other[,"V3.y"]) / sum(m.other[,"V3.x"]))

set.seed(25)
CNA.object <-CNA(genomdat = m.chr21[,"ratio"], chrom = m.chr21[,"V1"], maploc = m.chr21[,"V2"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, alpha = 0.0001, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2)

pdf(opt$output, height=4, width=7)
par(mar=c(2, 4, 1, 4))
plot(m.chr21$V2, m.chr21$ratio, ylim=c(-1.5,1.5), main=opt$patient, cex=0.8, pch=16, xlab="position chr21", ylab="log2 coverage ratio")
abline(0, 0)
for(i in 1:length(segs$output$loc.start)) {
	lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=6)
}
dev.off()
