# TABLE: sm_pathways.annotated

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("pathway file 1 not specified")
if (is.na(args[2])) stop("pathway file 2 not specified")

d <- read.csv(args[1], sep="\t", as.is=c(1, 10))
r <- read.csv(args[2], sep="\t", as.is=c(1, 10))

m <- merge(d, r, by=c("Pathway", "Name", "Class", "Size"), all=T, suffixes=c(".dia", ".rel"))
m <- m[order(m$p.value.rel),]

write.table(m, row.names=F, sep="\t")

