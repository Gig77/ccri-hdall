data <- read.csv("~/hdall/results/filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$status != "REJECT",]

patients <- levels(data$patient)
samples <- levels(data$sample)

min_cov = 30

#pdf("~/hdall/results/clonal-analysis/kernel-densities.allpatients.dia.pdf", height=12, paper='A4')
par(mfrow = c(6, 4), mar=c(2,0.5,2,0.5))
for(p in patients)
{
	d <- data$freq_leu[data$patient==p & data$sample=="rem_dia" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]
	plot(density(d), main=paste(p,"_dia (n=",length(d),")", sep=""), xlim=c(0,1), yaxt='n')
}
#dev.off()

pdf("~/hdall/results/clonal-analysis/kernel-densities.allpatients.rel.pdf", height=12, paper='A4')
par(mfrow = c(6, 4), mar=c(2,0.5,2,0.5))
for(p in patients)
{
	d <- data$freq_leu[data$patient==p & data$sample=="rem_rel" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]
	plot(density(d), main=paste(p,"_rel (n=",length(d),")", sep=""), xlim=c(0,1), yaxt='n')
}
dev.off()

pdf("~/hdall/results/clonal-analysis/patient715-kernel-density.pdf", width=12, paper='A4r')
par(mfrow = c(1, 2))
plot(density(data$freq_leu_norm[data$patient==715 & data$sample == "rem_dia" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]), main="715_dia")
plot(density(data$freq_leu_norm[data$patient==715 & data$sample == "rem_rel" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]), main="715_rel")
dev.off()

pdf("~/hdall/results/clonal-analysis/patient545-kernel-density.pdf", width=12, paper='A4r')
par(mfrow = c(1, 2))
plot(density(data$freq_leu_norm[data$patient==545 & data$sample == "rem_dia" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]), xlim=c(0,1), main="545_dia")
plot(density(data$freq_leu_norm[data$patient==545 & data$sample == "rem_rel" & t$dp_leu_tot >= min_cov & t$var_type == "snp"]), xlim=c(0,1), main="545_rel")
dev.off()
