options(warn=1)

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("ERROR: output file not specified")

var.dis <- read.delim("../filtered-variants.tsv", check.names=F, stringsAsFactor=F)
var.val <- read.delim("filtered-variants.reseq.tsv", check.names=F, stringsAsFactor=F)
panel <- read.delim("../panel-genes.tsv", check.names=F, stringsAsFactor=F, header=F)
panel[,2] <- "yes"
names(panel) <- c("gene", "panel")

var.dis <- var.dis[,1:37] # due to a bug in merge() function
var.val <- var.val[,1:37]
var.merged <- merge(var.dis, var.val, by=c("patient", "sample", "var_type", "chr", "pos", "dbSNP", "ref", "alt", "gene", "add_genes", "aa_change"), all=T, suffixes=c(".dis", ".val"))
var.merged <- merge(panel, var.merged, all.y=T)

write.table(var.merged, file=args[1], row.names=F, sep="\t", quote=F)
