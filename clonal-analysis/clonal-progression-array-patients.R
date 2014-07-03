#library("RColorBrewer")

min_cov = 50
patients <- c("C", "B", "818", "592", "446")  # 314, 399 excluded due to low coverage; patients 715, 545, 430 are main figure plots; 460 removed b/c only few mutations
#genes.to.label <- c("CREBBP", "KRAS", "NRAS", "FLT3", "PTPN11", "MLL2", "ATM", "WHSC1", "TRRAP", "TP53", "ZNF516", "NDC80", "TLX3", "GDPD2", "FBXL7")
genes.to.label <- c()
blast.count <- list(
		"C.dia" = 97, "C.rel" = 98, 
		"B.dia" = 94, "B.rel" = 90, 
		"460.dia" = 96, "460.rel" = 76, 
		"818.dia" = 95, "818.rel" = 89, 
		"592.dia" = 98, "592.rel" = 94, 
		"446.dia" = 85, "446.rel" = 94)

# exome seq variants
data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$status!="REJECT",]
data.filtered <- data[data$dp_leu_tot >= min_cov & data$dp_rem_tot >= min_cov & data$var_type == "snp" & data$'repeat'=="" & data$segdup=="",]
dia <- data.filtered[data.filtered$sample=="rem_dia", c("patient", "chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia)[7] <- "dia"
rel <- data.filtered[data.filtered$sample=="rem_rel", c("patient", "chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(rel)[7] <- "rel"

# reseq variants
data.reseq <- read.csv("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", sep="\t")
data.reseq <- data.reseq[data.reseq$patient == p & data.reseq$status!="REJECT",]
dia.reseq <- data.reseq[data.reseq$sample=="rem_dia", c("patient", "chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(dia.reseq)[7] <- "dia.reseq"
rel.reseq <- data.reseq[data.reseq$sample=="rem_rel", c("patient", "chr", "pos", "ref", "alt", "gene", "freq_leu_norm")]
names(rel.reseq)[7] <- "rel.reseq"

# merge both
m <- merge(dia, rel, all.x=T, all.y=T)
m <- merge(m, dia.reseq, all.x=T, all.y=T)
m <- merge(m, rel.reseq, all.x=T, all.y=T)

# use AF from resequencing data if available
m[!is.na(m$dia.reseq), "dia"] <- m[!is.na(m$dia.reseq), "dia.reseq"]
m[!is.na(m$rel.reseq), "rel"] <- m[!is.na(m$rel.reseq), "rel.reseq"]

# clean up AF (missing values, trimming to [0,1] interval)
m[is.na(m$dia), "dia"] <- 0
m[is.na(m$rel), "rel"] <- 0
m[m$dia>1, "dia"] <- 1
m[m$rel>1, "rel"] <- 1

# remove variants with AF < 10% at all timepoints (likely false positives)
m <- m[m$dia>=0.1 | m$rel>=0.1,]

# remove variants with AF > 70% (likely inaccurate measurement)
m <- m[m$dia<=0.7 & m$rel<=0.7,]

pdf("~/hdall/results/clonal-analysis/clonal-progression-array-patients.pdf")
par(mfrow = c(3, 2), mar=c(3,4,2,1))

for(p in patients)
{
	m.patient <- m[m$patient==p,]

	# assign colors based on kinetics
	m.patient$col <- "#000000"
	m.patient$col[m.patient$dia>0.15 & m.patient$rel>0.15] <- "#C6C7C8" # gray
	m.patient$col[m.patient$dia>0.25 & m.patient$rel==0] <- "#4A8ECC" # black
	m.patient$col[m.patient$dia==0 & m.patient$rel>0.25] <- "#BCAFD6" # violet
	m.patient$col[m.patient$dia==0 & m.patient$rel>0 & m.patient$rel<=0.25] <- "#FFF57C" # yellow
	m.patient$col[m.patient$dia > 0 & m.patient$dia<=0.25 & m.patient$rel > 0.25] <- "#F58E7D"
	m.patient$col[m.patient$dia > 0 & m.patient$dia<=0.25 & m.patient$rel==0] <- "#548A2F" # green
	
	plot(0, 0, xlim=c(1, 3.1), ylim=c(0, 0.7), type="n", xaxt="n", yaxt="n", xlab="", ylab="adj. allelic frequency", main=paste(p, " (n=", nrow(m.patient), ")", sep=""))
	axis(1, at=c(1.3, 2.8), labels=c(paste0("diagnosis\n(", blast.count[[paste0(p, ".dia")]], "% blasts)"), paste0("relapse\n(", blast.count[[paste0(p, ".rel")]], "% blasts)")), padj=0.5)
	axis(2, at = seq(0, 1, 0.1), las = 1) 
	
	for(i in 1:nrow(m.patient)) {
		fdia <- m.patient[i,"dia"]
		frel <- m.patient[i,"rel"]
		lines(c(1.3, 2.8), c(fdia, frel), type="b", col=m.patient[i,"col"], lwd=ifelse(m.patient[i, "gene"] %in% genes.to.label, 3, 1))
		
		if (m.patient[i, "gene"] %in% genes.to.label) {
			if (fdia > 0) { text(1.25, fdia, paste0(m.patient[i, "gene"], " -"), adj=c(1, 0.5), cex=0.8) }
			if (frel > 0) { text(2.85, frel, paste0("- ", m.patient[i, "gene"]), adj=c(0, 0.5), cex=0.8) }
		}
	}	
}

dev.off()