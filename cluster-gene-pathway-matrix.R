options(warn=1)

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("input file not specified")
if (is.na(args[2])) stop("output file not specified")

x <- as.matrix(read.csv(args[1], sep="\t", row.names=1, as.is=TRUE, check.names=F))
hc.rows <- hclust(dist(x, method="binary"), method="average")
hc.cols <- hclust(dist(t(x), method="binary"), method="average")
x <- x[hc.rows$order, hc.cols$order]
write.table(x, file=args[2], col.names=NA, row.names=T, sep="\t", quote=F)

