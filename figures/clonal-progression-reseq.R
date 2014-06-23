#library("RColorBrewer")
#X11.options(type="Xlib")

# 715, 545, Y, 592, 430
p <- "715"
min.cov <- 9
cov.max.std.dev <- 2
max.af <- 0.7
min.af <- 0.1
genes.to.label <- c("CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11", "MLL2")
exclude.chr <- c("chrX", "chrY")
blast.count <- list("715.dia" = 92, "715.rel" = 92, "715.rel3" = "92", "545.dia" = 97, "545.rel" = 86, "Y.dia" = 98, "Y.rel" = 95, "592.dia" = 94, "592.rel" = 98, "430.dia" = 81, "430.rel" = 72)

# exome seq variants
data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$patient == p & data$status!="REJECT",]
#data <- data[data$var_type=="snp",]  # variant type filter
data <- data[data$non_silent==1,]  # only non-silent
data <- data[!(data$chr %in% exclude.chr),]  # sex chromosome filter
data <- data[data$dp_leu_tot >= min.cov & data$dp_rem_tot >= min.cov,]  # minimum coverage filter
data <- data[data$dp_rem_tot < mean(data$dp_rem_tot) + cov.max.std.dev * sd(data$dp_rem_tot),] # maximum coverage filter

# patient 715, relapse 3
rel3 <- read.csv("/data/christian/p2ry8-crlf2/results/2014-03-27-old-reference/filtered-variants.cosmic.normaf.tsv", sep="\t", stringsAsFactor=F)
rel3 <- rel3[rel3$patient == p & rel3$status!="REJECT",]
#rel3 <- rel3[rel3$var_type=="snp",] # variant type filter
rel3 <- rel3[rel3$non_silent==1,] # only non-silent
rel3 <- rel3[!(rel3$chr %in% exclude.chr),] # sex chromosome filter
rel3 <- rel3[rel3$dp_leu_tot >= min.cov & rel3$dp_rem_tot >= min.cov,] # minimum coverage filter
rel3 <- rel3[rel3$dp_rem_tot < mean(rel3$dp_rem_tot) + cov.max.std.dev * sd(rel3$dp_rem_tot),] # maximum coverage filter

dia <- data[data$sample=="rem_dia", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia)[6] <- "dia"
rel <- data[data$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(rel)[6] <- "rel"
rel3 <- rel3[rel3$sample=="rem_rel3", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(rel3)[6] <- "rel3"

# reseq variants
data.reseq <- read.csv("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", sep="\t")
data.reseq <- data.reseq[data.reseq$patient == p & data.reseq$status!="REJECT",]
dia.reseq <- data.reseq[data.reseq$sample=="rem_dia", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia.reseq)[6] <- "dia.reseq"

if (p == "715") {
	rel.reseq <- data.reseq[data.reseq$sample=="rem_rel2", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
	rel3.reseq <- data.reseq[data.reseq$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
	names(rel3.reseq)[6] <- "rel3.reseq"
} else {
	rel.reseq <- data.reseq[data.reseq$sample=="rem_rel", c("chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
}
names(rel.reseq)[6] <- "rel.reseq"


# merge both
m <- merge(dia, rel, all.x=T, all.y=T)
m <- merge(m, rel3, all.x=T, all.y=T)
m <- merge(m, dia.reseq, all.x=T, all.y=T)
m <- merge(m, rel.reseq, all.x=T, all.y=T)

# use AF from resequencing data if available
# NOTE: patient Y has a much lower BLAST count at resequencing (~30%), because the sample was not taken from the BM but peripheral blood
if (p != "Y") {
	m[!is.na(m$dia.reseq), "dia"] <- m[!is.na(m$dia.reseq), "dia.reseq"]
	m[!is.na(m$rel.reseq), "rel"] <- m[!is.na(m$rel.reseq), "rel.reseq"]
}

if (p == "715") {
	m <- merge(m, rel3.reseq, all.x=T, all.y=T)
	m[!is.na(m$rel3.reseq), "rel3"] <- m[!is.na(m$rel3.reseq), "rel3.reseq"]
}

#par(mfrow = c(4, 5), mar=c(2,3,2,1))

m[is.na(m$dia), "dia"] <- 0
m[is.na(m$rel), "rel"] <- 0
m[is.na(m$rel3), "rel3"] <- 0
m[m$dia>1, "dia"] <- 1
m[m$rel>1, "rel"] <- 1
m[m$rel>1, "rel3"] <- 1

# remove variants with AF < 10% at all timepoints (likely false positives)
m <- m[m$dia>=min.af | m$rel>=min.af | m$rel3>=min.af,]

# remove variants with AF > 70% (likely inaccurate measurement)
m <- m[m$dia<=max.af & m$rel<=max.af & m$rel3<=max.af,]

if (p != "715")
{
	#png(paste0("~/hdall/results/figures/clonal-progression-", p, ".png"))
	pdf(paste0("~/hdall/results/figures/clonal-progression-", p, ".pdf"))
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
} else {
	m$col <- "#000000"
	m$col[m$dia>0.25 & m$rel>0.1 & m$rel3>0.25] <- "#E6E7E8"
	m$col[m$dia>0.25 & m$rel==0 & m$rel3==0] <- "#4A8ECC"
	m$col[m$dia>0.25 & m$rel>0 & m$rel<=0.25 & m$rel3==0] <- "#231F20"
	m$col[m$dia==0 & m$rel>0.25 & m$rel3 > 0.25] <- "#C9252B"
	m$col[m$dia==0 & m$rel>0.25 & m$rel3 == 0] <- "#76C58F"
	m$col[m$dia==0 & m$rel<=0.25 & m$rel3 == 0] <- "#548A2F"
	m$col[m$dia==0 & m$rel<=0.25 & m$rel3 > 0.25] <- "#F58E7D"
	m$col[m$dia==0 & m$rel==0 & m$rel3 > 0.25] <- "#BCAFD6"
	m$col[m$dia==0 & m$rel==0 & m$rel3 <= 0.25] <- "#FFF57C"
	
	#png(paste0("~/hdall/results/figures/clonal-progression-", p, ".png"), width=1024)
	pdf(paste0("~/hdall/results/figures/clonal-progression-", p, ".pdf"), width=12)
	plot(0, 0, xlim=c(1, 4.6), ylim=c(0, 0.7), type="n", xaxt="n", yaxt="n", xlab="", ylab="allelic frequency", main=paste(p, " (n=", nrow(m), ")", sep=""))
	axis(1, at=c(1.3, 2.8, 4.3), labels=c(paste0("diagnosis\n(", blast.count[[paste0(p, ".dia")]], "% blasts)"), paste0("relapse 1\n(", blast.count[[paste0(p, ".rel")]], "% blasts)"), paste0("relapse 3\n(", blast.count[[paste0(p, ".rel3")]], "% blasts)")), padj=0.5)
	axis(2, at = seq(0, 1, 0.1), las = 1) 
	
	for(i in 1:nrow(m)) {
		fdia <- m[i,"dia"]
		frel <- m[i,"rel"]
		frel3 <- m[i,"rel3"]
		lines(c(1.3, 2.8, 4.3), c(fdia, frel, frel3), type="b", col=m[i,"col"], lwd=ifelse(m[i, "gene"] %in% genes.to.label, 4, 1))
		
		if (m[i, "gene"] %in% genes.to.label) {
			text(1.25, fdia, paste0(m[i, "gene"], " --"), adj=c(1, 0.5), cex=0.8)
			text(4.35, frel3, paste0("-- ", m[i, "gene"]), adj=c(0, 0.5), cex=0.8)
		}
	}	
	
	write.table(m, file="~/hdall/results/figures/clonal-progression-715.data.tsv", col.names=T, row.names=F, sep="\t", quote=F)
}

dev.off()