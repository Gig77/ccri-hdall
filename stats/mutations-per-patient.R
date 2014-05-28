options(warn=1)
library(beeswarm)

# TABLE: filtered-variants.cosmic.normaf.tsv
sv <- read.delim("filtered-variants.cosmic.normaf.tsv")
sv <- subset(sv, patient != "E")

sv.pass <- subset(sv, status != "REJECT")
#sv.pass.nonsilent <- subset(sv.pass, non_silent == 1)
sv.pass.af10 <- subset(sv.pass, freq_leu >= 0.1)

t <- as.data.frame.matrix(table(sv.pass.af10$patient, sv.pass.af10$sample))

pdf("stats/mutations-per-patient-dia-vs-rel.pdf")
sig <- t.test(t$rem_dia, t$rem_rel, paired=T)
boxplot(t, col="cyan", outline=T, main="No. somatic mutations", names=c("diagnosis (n=19)", "relapse (n=19)"), xlab=paste("p=", format(sig$p.value, digits=1), ""))
beeswarm(t, pch=21, bg="grey", add=T)
dev.off()