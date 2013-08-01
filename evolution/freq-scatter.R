options(warn=1)

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("filtered variants not specified")
if (is.na(args[2])) stop("output PDF not specified")

# estimate minimal residual disease
patients <- c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "1021247", "A", "B", "C", "D", "E", "X", "Y") 
#patients <- c("399", "314")
t <- read.csv("~/hdall/results/filtered-variants.tsv", sep="\t")

pdf("~/hdall/results/evolution/freq-scatter.pdf", width=12, paper='A4r')
par(mfrow=c(4,5), mar=c(2,1,2,1))
for(p in patients) {
	freqd <- t[t$patient==p & t$sample=="rem_dia", c("chr", "pos", "freq")]
	colnames(freqd)[3] = "freq_dia"
	freqr <- t[t$patient==p & t$sample=="rem_rel", c("chr", "pos", "freq")]
	colnames(freqr)[3] = "freq_rel"
	merged <- merge(freqd, freqr, all.x=T, all.y =T)
	merged$freq_dia[is.na(merged$freq_dia)]	= 0
	merged$freq_rel[is.na(merged$freq_rel)]	= 0
	smoothScatter(merged$freq_dia, merged$freq_rel, nbin=100, bandwidth=0.01, xlim=c(0,1), ylim=c(0,1), main=paste(p, " (n=", length(merged$freq_dia), ")", sep=""))
	#plot(merged$freq_dia, merged$freq_rel, col="#00000033", xlim=c(0,1), ylim=c(0,1), main=paste(p, " (n=", length(merged$freq_dia), ")", sep=""))
}
dev.off()
