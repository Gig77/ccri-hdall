options(warn=1)
library(beeswarm)

# TABLE: filtered-variants.cosmic.normaf.tsv
sv <- read.delim("filtered-variants.cosmic.normaf.tsv")

sv.pass <- subset(sv, status != "REJECT")
sv.pass.nonsilent <- subset(sv.pass, effect %in% c("STOP_GAINED", "SPLICE_SITE_DONOR", "SPLICE_SITE_ACCEPTOR", "FRAME_SHIFT", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "NON_SYNONYMOUS_CODING", "CODON_INSERTION", "CODON_CHANGE_PLUS_CODON_DELETION", "NON_SYNONYMOUS_START", "START_LOST"))
sv.pass.nonsilent.af20 <- subset(sv.pass.nonsilent, freq_leu >= 0.2)

t <- as.data.frame.matrix(table(sv.pass.nonsilent.af20$patient, sv.pass.nonsilent.af20$sample))

pdf("stats/mutations-per-patient-dia-vs-rel.pdf")
sig <- t.test(t$rem_dia, t$rem_rel)
boxplot(t, col="cyan", outline=F, main="No. non-silent somatic mutations (AF>=20%)", names=c("diagnosis (n=20)", "relapse (n=20)"), xlab=paste("p=", format(sig$p.value, digits=1), ""))
beeswarm(t, pch=21, bg="grey", add=T)
dev.off()