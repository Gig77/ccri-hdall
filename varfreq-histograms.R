patients <- levels(data$patient)
samples <- levels(data$sample)

pdf("~/hdall/results/kernel-densities.allpatients.dia.pdf", height=12, paper='A4')
par(mfrow = c(6, 4), mar=c(2,0.5,2,0.5))
for(p in patients)
{
	d <- data$freq[data$patient==p & data$sample=="rem_dia"]
	plot(density(d), main=paste(p,"_dia (n=",length(d),")", sep=""), yaxt='n')
}
dev.off()

pdf("~/hdall/results/kernel-densities.allpatients.rel.pdf", height=12, paper='A4')
par(mfrow = c(6, 4), mar=c(2,0.5,2,0.5))
for(p in patients)
{
	d <- data$freq[data$patient==p & data$sample=="rem_rel"]
	plot(density(d), main=paste(p,"_rel (n=",length(d),")", sep=""), yaxt='n')
}
dev.off()

pdf("~/hdall/results/patient715-kernel-density.pdf", width=12, paper='A4r')
par(mfrow = c(1, 2))
plot(density(data$freq[data$patient==715 & data$sample == "rem_dia"]), main="715_dia")
plot(density(data$freq[data$patient==715 & data$sample == "rem_rel"]), main="715_rel")
dev.off()

pdf("~/hdall/results/patient545-kernel-density.pdf", width=12, paper='A4r')
par(mfrow = c(1, 2))
plot(density(data$freq[data$patient==545 & data$sample == "rem_dia"]), main="545_dia")
plot(density(data$freq[data$patient==545 & data$sample == "rem_rel"]), main="545_rel")
dev.off()
