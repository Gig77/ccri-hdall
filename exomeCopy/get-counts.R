library(exomeCopy)

target.file <- "~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr"
bam.files <- list.files(path="~/hdall/data/bam", pattern="*.bam$", full.names=T)
bam.files <- bam.files[!grepl("test", bam.files)]
sample.names <- paste0(sub("/home/STANNANET/christian.frech/hdall/data/bam/(.+).merged.duplicate_marked.realigned.recalibrated.bam", "\\1", bam.files))
reference.file <- "~/generic/data/hg19/ucsc.hg19.fasta"

target.df <- read.delim(target.file, header = FALSE)
target <- GRanges(seqname = target.df[, 1], IRanges(start = target.df[,2] + 1, end = target.df[, 3]))
counts <- RangedData(space = seqnames(target), ranges = ranges(target))

for (i in 1:length(bam.files)) {
	counts[[sample.names[i]]] <- countBamInGRanges(bam.files[i], target, read.width=100, remove.dup=FALSE)
}

names(counts@values@unlistData@listData) <- gsub("X(\\d)", "\\1", names(counts@values@unlistData@listData), perl=T)
save(counts, file="~/hdall/results/exomeCopy/counts.RData")

