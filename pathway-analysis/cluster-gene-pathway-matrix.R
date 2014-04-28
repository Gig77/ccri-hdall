options(warn=1)

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("input file not specified")
if (is.na(args[2])) stop("output file not specified")
#if (is.na(args[3])) stop("PDF output file for hierarchical gene tree not specified")
#if (is.na(args[4])) stop("PDF output file for hierarchical pathway tree not specified")

if (length(readLines(args[1])) == 0 || length(readLines(args[1])) == 1)
{
	print(paste("WARNING: Input file", args[1], "is empty. Cannot perform clustering of gene-pathway-matrix"))
	file.copy(args[1], args[2], overwrite=T)
	quit()
}
	 
x <- as.matrix(read.csv(args[1], sep="\t", row.names=1, as.is=TRUE, check.names=F))

if (ncol(x) < 2)
{
	print("WARNING: Input file args[1] contains less then two pathways. Cannot perform clustering of gene-pathway-matrix")
	write.table(x, file=args[2], col.names=NA, row.names=T, sep="\t", quote=F)
	quit()
}

# suppress leaf labels
#rownames(x) <- rep("", nrow(x))
#colnames(x) <- rep("", ncol(x))

hc.rows <- hclust(dist(x, method="binary"), method="average")
hc.cols <- hclust(dist(t(x), method="binary"), method="average")

x <- x[hc.rows$order, hc.cols$order]
write.table(x, file=args[2], col.names=NA, row.names=T, sep="\t", quote=F)

if (!is.na(args[3]) && nrow(x) > 2) {
	pdf(args[3], width=30, height=20)
	par(cex=0.6)
	plot(hc.rows, main=args[3], xlab=)
	dev.off()
}

if (!is.na(args[4]) && ncol(x) > 2) {
	pdf(args[4], width=30, height=20)
	par(cex=0.6)
	plot(hc.cols, main=args[4])
	dev.off()
}
