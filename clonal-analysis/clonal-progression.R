#library("RColorBrewer")

data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$status!="REJECT",]

#patients <- c("A")
patients <- levels(data$patient)
#palette <- colorRampPalette(brewer.pal(8, "Dark2"))(110)
min_cov = 30

pdf("~/hdall/results/clonal-analysis/clonal-progression.pdf", width=12)
par(mfrow = c(4, 5), mar=c(2,3,2,1))

for(p in patients)
{
	data.filtered <- data[data$patient==p & data$dp_leu_tot >= min_cov & data$dp_rem_tot >= min_cov & data$var_type == "snp" & data$'repeat'=="" & data$segdup=="",]
	dia <- data.filtered[data.filtered$sample=="rem_dia", c("chr", "pos", "ref", "alt", "freq_leu_norm")]
	names(dia)[5] <- "dia"
	rel <- data.filtered[data.filtered$sample=="rem_rel", c("chr", "pos", "ref", "alt", "freq_leu_norm")]
	names(rel)[5] <- "rel"
	m <- merge(dia, rel, all.x=T, all.y=T)
	m[is.na(m$dia), "dia"] <- 0
	m[is.na(m$rel), "rel"] <- 0
	m[m$dia>1, "dia"] <- 1
	m[m$rel>1, "rel"] <- 1
	
	plot(0, 0, xlim=c(1, 2), ylim=c(0, 1), type="n", xaxt="n", yaxt="n", xlab="", ylab="allelic frequency", main=paste(p, " (n=", nrow(m), ")", sep=""))
	axis(1, at=c(1, 2), labels=c("dia", "rel"))
	axis(2, at = seq(0, 1, 0.1), las = 1) 

	for(i in 1:nrow(m)) {
		fdia <- m[i,"dia"]
		frel <- m[i,"rel"]
		#lines(c(1, 2), c(min(m[i,"dia"], 1), min(m[i,"rel"], 1)), type="b", col=rgb(min(m[i,"dia"], 1),0,min(m[i,"rel"], 1), min(c(mean(c(m[i,"dia"], m[i,"rel"]))/2, 1))))
		lines(c(1, 2), c(fdia, frel), type="b", col=rgb(fdia, 0, frel, ((fdia+frel)/2)^2))
	}
}

dev.off()