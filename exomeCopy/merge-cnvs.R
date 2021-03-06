options(warn=1)
library(GenomicRanges)
library(reshape)

files <- list.files(path="/mnt/projects/hdall/results/exomeCopy", pattern="*.compiled-segments.tsv$", full.names=T)

# add HD ALL xenograft samples sequenced with the P2RY8-CRLF2 cohort
files <- c(files, list.files(path="/mnt/projects/p2ry8-crlf2/results/exomeCopy", pattern="m.*.compiled-segments.tsv$", full.names=T))

d <- read.delim(files[1], stringsAsFactor=F)
for (i in 2:length(files)) {
	d <- rbind(d, read.delim(files[i], stringsAsFactor=F))
}
d <- d[!d$sample.name %in% c("Y_rel", "X_rel", "m1957-715-rel"),] # ignore crappy samples
d$space <- as.factor(gsub("^(\\d+|X|Y)$", "chr\\1", d$space, perl=T))

gr <- GRanges(seqnames=d$space, ranges=IRanges(start=d$start, end=d$end), sample.name=d$sample.name, copy.count=d$copy.count, log.odds=d$log.odds, nranges=d$nranges, targeted.bp=d$targeted.bp, genes=d$genes)

# find overlaps
o <- findOverlaps(gr, gr)
o <- o[o@queryHits != o@subjectHits]
o <- o[(gr$copy.count[o@queryHits] > 2 & gr$copy.count[o@subjectHits] > 2) | (gr$copy.count[o@queryHits] < 2 & gr$copy.count[o@subjectHits] < 2)] 
o <- o[width(gr)[o@queryHits] >= 20 & width(gr)[o@subjectHits] >= 20] # remove overlaps with tiny single-exon segments

# determine overlap in percent of shared exons
ex <- read.delim("/mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr", header=F)
er <- GRanges(seqnames=ex$V1, ranges=IRanges(start=ex$V2, end=ex$V3))
oex <- findOverlaps(gr, er)
segex <- cast(as.data.frame(oex), formula=queryHits~., value="subjectHits", fun.aggregate=function(x) { paste(x, collapse=",") })
names(segex) <- c("subjectHits", "subjectExons")
om <- merge(as.data.frame(o), segex, all.x=T)
names(segex) <- c("queryHits", "queryExons")
om <- merge(om, segex, all.x=T)
om <- om[order(om$queryHits, om$subjectHits),]
o@metadata[["pct_shared_exons"]] <- apply(om, 1, function(x) { qex <- unlist(strsplit(x["queryExons"], ",", fixed=T)) ; sex <- unlist(strsplit(x["subjectExons"], ",", fixed=T)) ; length(intersect(qex, sex)) / length(union(qex, sex)) })
o <- o[o@metadata[["pct_shared_exons"]] >= 0.3]

gr$overlap.samples <- as.character(NA)
gr$overlap.count <- as.integer(NA)
gr$overlap.count.tumor <- as.integer(NA)
library(plyr)
o.aggregated <- ddply(as.data.frame(o), .(queryHits), function(x) { 
			c(overlap.samples=paste(unique(gr$sample.name[x$subjectHits]), collapse=","), 
			  overlap.count=length(unique(gr$sample.name[x$subjectHits])),
			  overlap.count.tumor=sum(!grepl("_rem$", unique(gr$sample.name[x$subjectHits])))) 
  })
gr$overlap.samples[o.aggregated$queryHits] <- o.aggregated$overlap.samples
gr$overlap.count[o.aggregated$queryHits] <- o.aggregated$overlap.count
gr$overlap.count.tumor[o.aggregated$queryHits] <- o.aggregated$overlap.count.tumor

# write table
write.table(as.data.frame(gr), file="/mnt/projects/hdall/results/exomeCopy/allpatients.compiled-segments.exomeCopy.tsv.part", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

