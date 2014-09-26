options(warn=1)
library(optparse)
library(gtools)

option_list <- list(
		make_option("--sample", type="character", help="sample name")
		)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- data.frame(sample="592_rel", stringsAsFactors=F)

if (invalid(opt$sample)) stop("sample not specified")

library(exomeCopy)

# ======================================
# combine individual chromosome to single virtual chromosome to improve optimization; 
# otherwise, CN estimates are off for single chromosomes with strong GC bias (e.g. chr13 in patient HV57R
#=======================================

combine_chr <- function(rd1, rd2) {
	print(paste0("Merging ", names(rd2), ": ", max(end(rd1))))
	tmp <- rd2
	ranges(tmp) <- shift(ranges(tmp), max(end(rd1)))
	names(tmp) <- names(rd1)
	rbind(rd1, tmp)
}

split_fit <- function(fit, chr, size) {
	extract <- end(fit@ranges[[1]]) <= size

	fit_first <- fit
	fit_first@ranges[[1]] <- fit_first@ranges[[1]][extract]
	fit_first@O.norm <- fit_first@O.norm[extract] 
	fit_first@log.odds <- fit_first@log.odds[extract] 
	fit_first@path <- fit_first@path[extract]
	names(fit_first@ranges) <- chr

	fit_second <- fit
	fit_second@ranges[[1]] <- fit_second@ranges[[1]][!extract]
	fit_second@ranges[[1]] <- shift(fit_second@ranges[[1]], -size) 
	fit_second@O.norm <- fit_second@O.norm[!extract] 
	fit_second@log.odds <- fit_second@log.odds[!extract] 
	fit_second@path <- fit_second@path[!extract]
	
	list(fit_first, fit_second)
}

# load counts and background estimates
load("counts.bg.RData")

# we do not fit individual chromosomes separately as suggested by vignette, but use larger genomic junks
# gives better results, probably due to more accurate GC bias estimates
# however, we cannot do the whole genome at once, because of integer overflow (max=2^31-1) :-(
counts.combined.1 <- counts["chr1"]
for (chr in c("chr3", "chr5", "chr7", "chr9", "chr11", "chr13", "chr15", "chr17", "chr19", "chr21", "chrY")) {
	counts.combined.1 <- combine_chr(counts.combined.1, counts[chr])
}
counts.combined.2 <- counts["chr2"]
for (chr in c("chr4", "chr6", "chr8", "chr10", "chr12", "chr14", "chr16", "chr18", "chr20", "chr22", "chrX")) {
	counts.combined.2 <- combine_chr(counts.combined.2, counts[chr])
}

fit.all.1 <- exomeCopy(counts.combined.1["chr1"], opt$sample, X.names = c("log.bg", "GC", "GC.sq", "width"), S = 0:6, d = 2)
fit.all.2 <- exomeCopy(counts.combined.2["chr2"], opt$sample, X.names = c("log.bg", "GC", "GC.sq", "width"), S = 0:6, d = 2)

# split into chromosomes again
fit <- list()
fit[[opt$sample]] <- list()
fit.rest <- fit.all.1
for (chr in c("chr1", "chr3", "chr5", "chr7", "chr9", "chr11", "chr13", "chr15", "chr17", "chr19", "chr21", "chrY")) {
	splitted <- split_fit(fit.rest, chr, max(end(counts[chr])))
	fit[[opt$sample]][[chr]] <- splitted[[1]]
	fit.rest <- splitted[[2]] 
}
fit.rest <- fit.all.2
for (chr in c("chr2", "chr4", "chr6", "chr8", "chr10", "chr12", "chr14", "chr16", "chr18", "chr20", "chr22", "chrX")) {
	splitted <- split_fit(fit.rest, chr, max(end(counts[chr])))
	fit[[opt$sample]][[chr]] <- splitted[[1]]
	fit.rest <- splitted[[2]] 
}

# plot results
pdf(paste0(opt$sample, ".combined.pdf"), width=15, height=7)
for (chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
	plot(fit[[opt$sample]][[chr]], main=chr, cex=0.2, ylim=c(0, 6), xlab='', ylab='', cex.main=2, cex.axis=1, col=c("red", "orange", "gray", "deepskyblue", "blue", "blue2", "blue4")) 
}
dev.off()

# get compiled segments
compiled.segments <- compileCopyCountSegments(fit)
compiled.segments$genes <- as.character(NA)

# annotate segments with overlapping gene names
if(file.exists("genes.GRCh37v75.biomart.RData")) {
	load("genes.GRCh37v75.biomart.RData")
} else {
	library("biomaRt")
	mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75
	genes <- getBM(c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"), mart=mart)
	save(genes, file="genes.GRCh37v75.biomart.RData")
}
gr.genes <- GRanges(seqnames=paste0("chr", genes$chromosome_name), ranges=IRanges(start=genes$start_position, end=genes$end_position), names=genes$hgnc_symbol)

o <- findOverlaps(compiled.segments, gr.genes)
for(i in 1:nrow(compiled.segments)) { gnames <- values(gr.genes[o@subjectHits[o@queryHits==i]])$names; compiled.segments[["genes"]][i] <- paste(gnames[gnames!=""], collapse=",") }

# write table
write.table(compiled.segments, file=paste0(opt$sample, ".compiled-segments.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

