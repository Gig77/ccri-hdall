p <- read.delim("/mnt/projects/hdall/results/filtered-variants.cosmic.tsv")
x <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.xenografts.tsv")
p$chr <- gsub("chr", "", as.character(p$chr))
p$chr[p$chr=="M"] <- "MT"
p.715r1 <- p[p$patient=="715" & p$sample=="rem_rel" & p$status=="PASS",c("chr","pos","ref","alt","gene","effect","aa_change","freq_leu")]
x.715r1 <- x[x$patient=="715" & x$sample=="rem_rel_xeno_m1957" & x$status=="PASS",c("chr","pos","ref","alt","gene","effect","freq_leu")]
m <- merge(p.715r1, x.715r1, by=c("chr", "pos", "ref", "alt", "gene", "effect"), suffixes=c(".primary", ".xeno"), all=T)
write.table(m, file="/mnt/projects/hdall/results/check-xeno-af.tsv", quote=F, sep="\t", row.names=F)
