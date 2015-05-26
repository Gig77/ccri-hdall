options(warn=1)
#library(vioplot)

gs <- read.delim("/mnt/projects/hdall/results/gene-size.txt", stringsAsFactors=F)
mm <- read.delim("/mnt/projects/hdall/results/gene-patient-matrix.annotated.tsv", check.names=F, stringsAsFactors=F)
il <- read.delim("/mnt/projects/generic/data/illumina/truseq_exome_targeted_regions.hg19.bed.chr", check.names=F, stringsAsFactors=F)

il.genes <- data.frame(HGNC=unique(matrix(unlist(strsplit(il[,4], ":", fixed=T)), ncol=3, byrow=T)[,3]), kit="yes", stringsAsFactors=F)
names(mm)[1] <- "HGNC"

gs.kit <- merge(gs[gs$cds_len!="NA",], il.genes)

pdf("/mnt/projects/hdall/results/stats/gene-length-bias.pdf", width=10)
par(mfrow=c(1, 2), oma=c(1.5,0,1.5,0))

# check for significant difference between mutated and not-mutated genes
merged <- merge(gs.kit, mm[,c("HGNC", "tot-dia", "tot-rel", "tot-cons")], all.x=T)
merged[is.na(merged)] <- 0
#num.mut <- merged$'tot-dia'+merged$'tot-rel'-merged$'tot-cons'
num.mut.factor <- num.mut <- merged$'tot-rel'
num.mut.factor[num.mut.factor>0] = "mutated"
num.mut.factor[num.mut.factor!="mutated"] = "not mutated"
num.mut.factor <- factor(num.mut.factor, c("not mutated", "mutated"))
sig <- t.test(merged[num.mut==0, "cds_len"], merged[num.mut>0, "cds_len"])
plot(merged$cds_len~num.mut.factor, ylim=c(0,6000), outline=F, xlab=paste("p=", format(sig$p.value, digits=1), sep=""), ylab="CDS length (bp)", col="cyan")
text(1:2, rep(5000,2), c(paste("n=", sum(num.mut==0), sep=""), paste("n=", sum(num.mut>0), sep="")))
#vioplot(merged$cds_len[num.mut.factor=="0"], merged$cds_len[num.mut.factor=="1"], merged$cds_len[num.mut.factor=="2+"], 
#		ylim=c(0, 20000), 
#		names=c(paste("0 (n=", sum(num.mut==0), ")", sep=""), paste("1 (n=", sum(num.mut==1), ")", sep=""), paste("2+ (n=", sum(num.mut>=2), ")", sep="")),
#		col="gray")
#mtext("number somatic mutations", 1, line=3)

# check for significant difference between MuSiC significant genes and non-significant genes
#mm.rel <- mm[mm$'tot-rel' > 0 & !is.na(mm$cds_len), c("cds_len", "p-gene-rel")]
merged <- merge(gs.kit, mm[, c("HGNC", "p-gene-rel")], all.x=T)
num.mut <- merged$'p-gene-rel'
num.mut[is.na(num.mut)] <- 1
num.mut[num.mut<0.05] <- "MuSiC significant"
num.mut[num.mut!="MuSiC significant"] <- "MuSiC not significant"
num.mut <- as.factor(num.mut)
sig <- t.test(merged[num.mut=="MuSiC significant", "cds_len"], merged[num.mut=="MuSiC not significant", "cds_len"])
plot(merged$cds_len~num.mut, ylim=c(0,6000), outline=F, xlab=paste("p=", format(sig$p.value, digits=1), sep=""), ylab="CDS length (bp)", col="cyan")
text(1:2, rep(5000,2), c(paste("n=", sum(num.mut=="MuSiC not significant"), sep=""), paste("n=", sum(num.mut=="MuSiC significant"), sep="")))
#points(rep(1, sum(num.mut=="MuSiC not significant")), mm.rel$cds_len[num.mut=="MuSiC not significant"], col=rgb(0,0,0,0.05))
#points(rep(2, sum(num.mut=="MuSiC significant")), mm.rel$cds_len[num.mut=="MuSiC significant"], col=rgb(0,0,0,0.05))

mtext("Average gene size by mutation status", side=3, outer=T, cex=1.5, line=-1)

dev.off()