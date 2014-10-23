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

# split by sex
sex <- read.delim("~/hdall/results/patient_sex.tsv")
sample.names.bg.male <- sample.names.bg[sample.names.bg %in% paste0(sex$patient[sex$sex=="m"], "_rem")] 
sample.names.bg.female <- sample.names.bg[sample.names.bg %in% paste0(sex$patient[sex$sex=="f"], "_rem")] 

# exclude problematic samples
sample.names.bg <- sample.names.bg[!sample.names.bg %in% c("399_rem", "314_rem")]

counts[["GC"]] <- getGCcontent(target, reference.file)
counts[["GC.sq"]] <- counts$GC^2
counts[["bg"]] <- generateBackground(sample.names.bg, counts, median)
counts[["log.bg"]] <- log(counts[["bg"]] + 0.1)
counts[["bg.var"]] <- generateBackground(sample.names.bg, counts, var) 
counts[["bg.male"]] <- generateBackground(sample.names.bg.male, counts, median)
counts[["log.bg.male"]] <- log(counts[["bg.male"]] + 0.1)
counts[["bg.var.male"]] <- generateBackground(sample.names.bg.male, counts, var) 
counts[["log.bg.var.male"]] <- log(counts[["bg.var.male"]] + 0.1) 
counts[["bg.female"]] <- generateBackground(sample.names.bg.female, counts, median)
counts[["log.bg.female"]] <- log(counts[["bg.female"]] + 0.1)
counts[["bg.var.female"]] <- generateBackground(sample.names.bg.female, counts, var) 
counts[["log.bg.var.female"]] <- log(counts[["bg.var.female"]] + 0.1) 
counts[["width"]] <- width(counts)

save(counts, file="counts.bg.RData")

