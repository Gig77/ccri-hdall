library(exomeCopy)

target.file <- "~/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr"
bam.files <- list.files(path="~/hdall/data/bam", pattern="*.bam$", full.names=T)
bam.files <- bam.files[!grepl("test", bam.files)]
sample.names <- paste0(sub("/home/STANNANET/christian.frech/hdall/data/bam/(.+).merged.duplicate_marked.realigned.recalibrated.bam", "\\1", bam.files))
reference.file <- "~/generic/data/hg19/ucsc.hg19.fasta"

target.df <- read.delim(target.file, header = FALSE)
target <- GRanges(seqname = target.df[, 1], IRanges(start = target.df[,2] + 1, end = target.df[, 3]))

# load counts
load("counts.RData")

# only remission samples
sample.names.bg <- sample.names[grep("rem$", sample.names)]

# exclude problematic samples
sample.names.bg <- sample.names.bg[!sample.names.bg %in% c("430_rem", "460_rem", "399_rem", "314_rem")]

counts[["GC"]] <- getGCcontent(target, reference.file)
counts[["GC.sq"]] <- counts$GC^2
counts[["bg"]] <- generateBackground(sample.names.bg, counts, median)
counts[["log.bg"]] <- log(counts[["bg"]] + 0.1)
counts[["bg.var"]] <- generateBackground(sample.names.bg, counts, var) 
counts[["width"]] <- width(counts)

# remove outliers; gives more robust parameter optimization results
counts <- counts[counts[["bg.var"]]>quantile(counts[["bg.var"]], 0.01) & counts[["bg.var"]]<quantile(counts[["bg.var"]], 0.99),]

save(counts, file="counts.bg.RData")

