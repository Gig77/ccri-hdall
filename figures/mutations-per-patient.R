library("RColorBrewer")
warnings()

patients <- c("1021247", "446", "818", "A", "592", "430", "792", "786", "B", "460", "545", "715", "X", "Y", "842", "C", "D", "399", "314")
effect.cds <- c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "START_GAINED", "START_LOST", "STOP_GAINED", "SYNONYMOUS_CODING", "SYNONYMOUS_STOP")

#patients <- c("314", "399") 

# TABLE: filtered-variants.tsv
# read and filter input data
t <- read.delim("~/hdall/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$freq_leu >= 0.1 & t$non_silent==1,]

# join diagnosis and relapse variants
t.dia <- t[t$sample == "rem_dia", c("patient", "chr", "pos", "ref", "alt", "freq_leu")]
names(t.dia)[6] <- "dia"
t.rel <- t[t$sample == "rem_rel", c("patient", "chr", "pos", "ref", "alt", "freq_leu")]
names(t.rel)[6] <- "rel"
m <- merge(t.dia, t.rel, all.x=T, all.y=T)

bardata <- matrix(rep(0, 3*length(patients)), nrow=length(patients))
colnames(bardata) <- c("relapse-specific", "diagnosis-specific", "conserved")
rownames(bardata) <- patients

# count mutations per patient
for(p in patients) {
	m.p <- m[m$patient==p,]
	bardata[p,1] <- sum(is.na(m.p$dia) & !is.na(m.p$rel))
	bardata[p,2] <- sum(!is.na(m.p$dia) & is.na(m.p$rel))
	bardata[p,3] <- sum(!is.na(m.p$dia) & !is.na(m.p$rel))
}
bardata <- bardata[order(rowSums(bardata)),]

make_plot <- function() {
	barplot(t(bardata), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, main="Non-silent mutations per patient", ylab="# non-silent somatic mutations (allelic frequency >= 10%)", ylim=c(0, 100))
	abline(h=seq(10, 100, by=10), col="gray")
	box()
	barplot(t(bardata), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, main="Non-silent mutations per patient", ylab="# non-silent somatic mutations (allelic frequency >= 10%)", ylim=c(0, 100), add=T)
	legend(x=1, y=95, bg="white", rev(colnames(bardata)), fill=rev(brewer.pal(3, "Set2")), cex=1.5)
}

# plot PNG
png("~/hdall/results/figures/mutation-per-patient.png", width=1600, height=1600, pointsize=40)
make_plot()
dev.off()

# plot PDF
pdf("~/hdall/results/figures/mutation-per-patient.pdf")
make_plot()
dev.off()

# box plot, statistical test
png("~/hdall/results/figures/nonsilent-mutations-dia-vs-rel.png", width=800, height=800, res=100)
mut <- as.data.frame(bardata)
mut$dia <- mut$'diagnosis-specific' + mut$conserved
mut$rel <- mut$'relapse-specific' + mut$conserved
test <- kruskal.test(list(mut$dia, mut$rel))
boxplot(mut$dia, mut$rel, xlab="timepoint", ylab="number non-silent mutations", na.action=na.exclude, outline=F, cex.axis=0.8, main=paste0("p=", sprintf("%.2g", test$p.value), " (Kruskal-Wallis)"), xaxt="n", cex.lab=1.3)
axis(1, at=c(1,2), labels=c("diagnosis", "relapse"))
stripchart(list(mut$dia, mut$rel), method="jitter", vertical=T, pch=19, col=c("red", "blue"), add=T)
dev.off()
