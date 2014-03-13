warnings()

library(vcd)
library(ggplot2)

min.af <- 0.10
min.dp.leu <- 10
min.dp.rem <- 10

# diagnosis, relapse 1
kamilla <- read.delim("~/hdall/results/filtered-variants.cosmic.normaf.tsv", stringsAsFactors=F, na.strings=c("NA", "n/d"))
maria <- read.delim("~/p2ry8-crlf2/results/filtered-variants.cosmic.normaf.tsv", stringsAsFactor=F, na.strings=c("NA", "n/d"))

m <- kamilla[kamilla$patient=="715" & kamilla$status != "REJECT" & kamilla$freq_leu >= min.af,]
m <- rbind(m, maria[maria$patient=="715" & maria$status != "REJECT" & maria$freq_leu >= min.af,])

pdf("~/hdall/results/patient715/patient715.pdf")

# number of variants per sample
barplot(table(m$sample), ylim=c(0,200))

# difference snps / indels
table.indel <- table(m[,c("sample", "var_type")])
p <- fisher.test(table.indel)$p.value
mosaic(table.indel, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.indel)(table.indel)

# difference silent vs. non-silent mutations in samples?
table.silent.nonsilent <- table(m[,c("sample", "non_silent")])
p <- fisher.test(table.silent.nonsilent)$p.value
mosaic(table.silent.nonsilent, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.silent.nonsilent)(table.silent.nonsilent)

# difference in effect?
table.effect <- table(m[m$effect %in% c("SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING"),c("sample", "effect")])
p <- fisher.test(table.effect)$p.value
mosaic(table.effect, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.effect)(table.effect)

# difference in impact?
table.impact <- table(m[,c("sample", "impact")])
p <- fisher.test(table.impact)$p.value
mosaic(table.impact, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.impact)(table.impact)

# difference in location?
m$gene.context <- NA
m$gene.context[m$effect=="INTRON"] <- "intronic"
m$gene.context[m$effect %in% c("INTERGENIC", "DOWNSTREAM", "UPSTREAM")] <- "intergenic" 
m$gene.context[m$effect %in% c("UTR_3_PRIME", "SPLICE_SITE_DONOR", "NON_SYNONYMOUS_CODING", "UTR_5_PRIME", "SYNONYMOUS_CODING", "START_GAINED", "FRAME_SHIFT", "STOP_GAINED", "EXON")] <- "genic"

table.gene.context <- table(m[,c("sample", "gene.context")])
p <- fisher.test(table.gene.context)$p.value
mosaic(table.gene.context, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.gene.context)(table.gene.context)

# difference in deleteriousness?
table.deleterious <- table(m[,c("sample", "deleterious")])
p <- fisher.test(table.deleterious)$p.value
mosaic(table.deleterious, pop=F, main=sprintf("p=%.2g", p))
labeling_cells(text=table.deleterious)(table.deleterious)

# difference in allelic frequency
p <- kruskal.test(sample~freq_leu, data=m, na.action=na.exclude)$p.value
boxplot(freq_leu~sample, data=m, main=sprintf("p=%.2g", p), ylab="allelic frequency", na.action=na.exclude, outline=F)
stripchart(freq_leu~sample, data=m, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=1:length(levels(as.factor(m$sample))), add=T)

# kernel density distribution of allelic frequencies by sample
ggplot(m, aes(freq_leu, fill = sample)) + geom_density(alpha = 0.2) + scale_x_continuous(limits = c(-0.2,1.2), breaks=seq(0, 1, 0.1))

# kernel density distribution of allelic frequencies by sample (normalized by chromosome number)
ggplot(m, aes(freq_leu_norm, fill = sample)) + geom_density(alpha = 0.2) + scale_x_continuous(limits = c(-0.2,1.2), breaks=seq(0, 1, 0.1))

# mutations per chromosome
table.chrom <- table(m[!m$chr %in% c("chr6_ssto_hap7", "chrUn_gl000220"),c("sample", "chr")])
p <- fisher.test(table.chrom, simulate.p.value=T, B=1e5)$p.value
table.chrom[1,] <- table.chrom[1,] / sum(table.chrom[1,])
table.chrom[2,] <- table.chrom[2,] / sum(table.chrom[2,])
table.chrom[3,] <- table.chrom[3,] / sum(table.chrom[3,])
barplot(table.chrom, beside=T, main=sprintf("p=%.2g", p), ylab="relative frequency", las=3)
legend(x=50, y=0.11, rownames(table.chrom), fill=gray.colors(nrow(table.chrom)))

# mutations per ploidy
table.ploidy <- table(m[,c("sample", "copy_no")])
p <- fisher.test(table.ploidy)$p.value
barplot(table.ploidy, beside=T, main=sprintf("p=%.2g", p), ylab="number mutations", xlab="ploidy")
legend(x=10, y=150, rownames(table.ploidy), fill=gray.colors(nrow(table.ploidy)))

# transitions/transversions all three samples
m$mut_type <- NA
m$mut_type[m$ref=="A" & m$alt=="G" | m$ref=="G" & m$alt=="A"] <- "A->G" 
m$mut_type[m$ref=="C" & m$alt=="T" | m$ref=="T" & m$alt=="C"] <- "C->T" 
m$mut_type[m$ref=="G" & m$alt=="T" | m$ref=="T" & m$alt=="G"] <- "G->T" 
m$mut_type[m$ref=="A" & m$alt=="C" | m$ref=="C" & m$alt=="A"] <- "A->C" 
m$mut_type[m$ref=="G" & m$alt=="C" | m$ref=="C" & m$alt=="G"] <- "C->G" 
m$mut_type[m$ref=="A" & m$alt=="T" | m$ref=="T" & m$alt=="A"] <- "A->T" 
table.muttype <- table(m[,c("sample", "mut_type")])
p <- fisher.test(table.muttype, simulate.p.value=T, B=1e5)$p.value

table.muttype[1,] <- table.muttype[1,] / sum(table.muttype[1,])
table.muttype[2,] <- table.muttype[2,] / sum(table.muttype[2,])
table.muttype[3,] <- table.muttype[3,] / sum(table.muttype[3,])

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(table.muttype), col=brewer.pal(6, "Set1"), xlab="sample", ylab="relative frequency", main=sprintf("p=%.2g", p))
legend("topright", inset=c(-0.2,0), legend=colnames(table.muttype), fill=rev(brewer.pal(6, "Set1")))

# transitions/transversions all three samples
table.muttype <- table(m[m$sample %in% c("rem_rel", "rem_rel3"), c("sample", "mut_type")])
p <- fisher.test(table.muttype, simulate.p.value=T, B=1e5)$p.value

table.muttype[1,] <- table.muttype[1,] / sum(table.muttype[1,])
table.muttype[2,] <- table.muttype[2,] / sum(table.muttype[2,])

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
barplot(t(table.muttype), col=brewer.pal(6, "Set1"), xlab="sample", ylab="relative frequency", main=sprintf("p=%.2g", p))
legend("topright", inset=c(-0.2,0), legend=colnames(table.muttype), fill=rev(brewer.pal(6, "Set1")))

dev.off()