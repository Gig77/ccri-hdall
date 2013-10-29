data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")

patients <- c("430")
#patients <- levels(data$patient)
samples <- levels(data$sample)

min_cov = 30

#par(mfrow = c(6, 4), mar=c(2,0.5,2,0.5))
for(p in patients)
{
	d <- data[data$patient==p & data$dp_leu_tot >= min_cov & data$var_type == "snp", c("chr", "pos", "sample", "freq_leu")]
	dia <- data[data$patient==p & data$sample=="rem_dia" & data$dp_leu_tot >= min_cov & data$var_type == "snp", c("chr", "pos", "freq_leu")]
	names(dia)[3] <- "dia"
	rel <- data[data$patient==p & data$sample=="rem_rel" & data$dp_leu_tot >= min_cov & data$var_type == "snp", c("chr", "pos", "freq_leu")]
	names(rel)[3] <- "rel"
	m <- merge(dia, rel, all.x=T, all.y=T)
	m[is.na(m$dia), "dia"] <- 0
	m[is.na(m$rel), "rel"] <- 0
	plot(d, type="b")
	#plot(density(d), main=paste(p,"_dia (n=",length(d),")", sep=""), xlim=c(0,1), yaxt='n')
}
