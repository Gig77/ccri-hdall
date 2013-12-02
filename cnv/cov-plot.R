dia <- read.delim("~/hdall/results/cnv/592_dia.coverage.tsv", check.names=F, stringsAsFactor=F, header=F)
dia.mito <- read.delim("~/hdall/results/cnv/592_dia.coverage.mito.tsv", check.names=F, stringsAsFactor=F, header=F)
rel <- read.delim("~/hdall/results/cnv/592_rel.coverage.tsv", check.names=F, stringsAsFactor=F, header=F)
rel.mito <- read.delim("~/hdall/results/cnv/592_rel.coverage.mito.tsv", check.names=F, stringsAsFactor=F, header=F)
rem <- read.delim("~/hdall/results/cnv/592_rem.coverage.tsv", check.names=F, stringsAsFactor=F, header=F)
rem.mito <- read.delim("~/hdall/results/cnv/592_rem.coverage.mito.tsv", check.names=F, stringsAsFactor=F, header=F)

# normalization factors
nf.dia = sum(dia.mito$V3)/sum(rem.mito$V3) 
nf.rel = sum(rel.mito$V3)/sum(rem.mito$V3) 

names(dia) <- c("chr", "start", "end", "name", "length", "strand", "total.dia", "avg.dia")
names(rem) <- c("chr", "start", "end", "name", "length", "strand", "total.rem", "avg.rem")
names(rel) <- c("chr", "start", "end", "name", "length", "strand", "total.rel", "avg.rel")

m <- merge(rem, dia, by=c("chr", "start", "end", "name", "length", "strand"))
m <- merge(m, rel, by=c("chr", "start", "end", "name", "length", "strand"))

# normalize coverage to total coverage
mn <- m
mn <- mn[mn$avg.dia >= 30 | mn$avg.rel >= 30 | mn$avg.rem >= 30,]
mn$avg.dia <- mn$avg.dia / sum(mn[mn$chr=="chr19", "avg.dia"]) * sum(mn[mn$chr=="chr19","avg.rem"])
mn$avg.rel <- mn$avg.rel / sum(mn[mn$chr=="chr19", "avg.rel"]) * sum(mn[mn$chr=="chr19","avg.rem"])
#mn$avg.dia <- mn$avg.dia / nf.dia
#mn$avg.rel <- mn$avg.rel / nf.rel

chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

pdf("~/hdall/results/cnv/cov-plot.592.pdf")
for (c in chr)
{
	par(mfrow=c(2,1))
	plot(mn[mn$chr==c,"start"],log2(mn[mn$chr==c,"avg.dia"]/mn[mn$chr==c,"avg.rem"]), ylim=c(-2, 2), col=rgb(0,0,0,0.05), main=paste("592", c, "dia"), ylab="log fold cov", xlab="", cex=0.2)
	plot(mn[mn$chr==c,"start"],log2(mn[mn$chr==c,"avg.rel"]/mn[mn$chr==c,"avg.rem"]), ylim=c(-2, 2), col=rgb(0,0,0,0.05), main=paste("592", c, "rel"), ylab="log fold cov", xlab="", cex=0.2)
}
dev.off()

