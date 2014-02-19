#library("RColorBrewer")

# 715, 545, Y, 592
p <- "592"
min.cov <- 30
cov.max.std.dev <- 2
genes.to.label <- c("CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11")
exclude.chr <- c("chrX", "chrY")
blast.count <- list("715.dia" = 92, "715.rel" = 92, "545.dia" = 97, "545.rel" = 86, "Y.dia" = 98, "Y.rel" = 95, "592.dia" = 94, "592.rel" = 98)

# exome seq variants
data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$patient == p & data$status!="REJECT" & data$dp_leu_tot >= min.cov & data$dp_rem_tot >= min.cov & !(data$chr %in% exclude.chr),]
data <- data[data$dp_rem_tot < mean(data$dp_rem_tot) + cov.max.std.dev * sd(data$dp_rem_tot),]
dia <- data[data$sample=="rem_dia", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia)[6] <- "dia"
rel <- data[data$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(rel)[6] <- "rel"

# reseq variants
data.reseq <- read.csv("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", sep="\t")
data.reseq <- data.reseq[data.reseq$patient == p & data.reseq$status!="REJECT",]
dia.reseq <- data.reseq[data.reseq$sample=="rem_dia", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia.reseq)[6] <- "dia.reseq"
if (p == "715") {
	rel.reseq <- data.reseq[data.reseq$sample=="rem_rel2", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
} else {
	rel.reseq <- data.reseq[data.reseq$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
}
names(rel.reseq)[6] <- "rel.reseq"

# merge both
m <- merge(dia, rel, all.x=T, all.y=T)
m <- merge(m, dia.reseq, all.x=T, all.y=T)
m <- merge(m, rel.reseq, all.x=T, all.y=T)

# patient Y has a much lower BLAST count at resequencing (~30%), because the sample was not taken from the BM but peripheral blood
if (p != "Y") {
	m[!is.na(m$dia.reseq), "dia"] <- m[!is.na(m$dia.reseq), "dia.reseq"]
	m[!is.na(m$rel.reseq), "rel"] <- m[!is.na(m$rel.reseq), "rel.reseq"]
}

png(paste0("~/hdall/results/figures/clonal-progression-", p, ".png"))
#par(mfrow = c(4, 5), mar=c(2,3,2,1))

m[is.na(m$dia), "dia"] <- 0
m[is.na(m$rel), "rel"] <- 0
m[m$dia>1, "dia"] <- 1
m[m$rel>1, "rel"] <- 1
	
plot(0, 0, xlim=c(1, 3.1), ylim=c(0, 0.7), type="n", xaxt="n", yaxt="n", xlab="", ylab="allelic frequency", main=paste(p, " (n=", nrow(m), ")", sep=""))
axis(1, at=c(1.3, 2.8), labels=c(paste0("diagnosis\n(", blast.count[[paste0(p, ".dia")]], "% blasts)"), paste0("relapse\n(", blast.count[[paste0(p, ".rel")]], "% blasts)")), padj=0.5)
axis(2, at = seq(0, 1, 0.1), las = 1) 

for(i in 1:nrow(m)) {
	fdia <- m[i,"dia"]
	frel <- m[i,"rel"]
	#lines(c(1, 2), c(min(m[i,"dia"], 1), min(m[i,"rel"], 1)), type="b", col=rgb(min(m[i,"dia"], 1),0,min(m[i,"rel"], 1), min(c(mean(c(m[i,"dia"], m[i,"rel"]))/2, 1))))
	lines(c(1.3, 2.8), c(fdia, frel), type="b", col=rgb(fdia, 0, frel, ifelse(m[i, "gene"] %in% genes.to.label, 1, 0.3)), lwd=ifelse(m[i, "gene"] %in% genes.to.label, 4, 1))
	
	if (m[i, "gene"] %in% genes.to.label) {
		if (fdia > 0) { text(1.25, fdia, paste0(m[i, "gene"], " --"), adj=c(1, 0.5), cex=0.8) }
		if (frel > 0) { text(2.85, frel, paste0("-- ", m[i, "gene"]), adj=c(0, 0.5), cex=0.8) }
	}
}

dev.off()