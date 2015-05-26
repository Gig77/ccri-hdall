library("RColorBrewer")
warnings()

patients <- c("1021247", "446", "818", "A", "592", "430", "792", "786", "B", "460", "545", "715", "X", "Y", "842", "C", "D", "399", "314")

# TABLE: filtered-variants.tsv
# read and filter input data
t <- read.delim("/mnt/projects/hdall/results/filtered-variants.tsv", stringsAsFactors=F)
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

# bar data
for(p in patients) {
	m.p <- m[m$patient==p,]
	bardata[p,1] <- sum(is.na(m.p$dia) & !is.na(m.p$rel))
	bardata[p,2] <- sum(!is.na(m.p$dia) & is.na(m.p$rel))
	bardata[p,3] <- sum(!is.na(m.p$dia) & !is.na(m.p$rel))
}
bardata <- bardata[order(rowSums(bardata)),]

# boxplot data
mut <- as.data.frame(bardata)
mut$dia <- mut$'diagnosis-specific' + mut$conserved
mut$rel <- mut$'relapse-specific' + mut$conserved
test <- kruskal.test(list(mut$dia, mut$rel))

# plot
make_plot <- function() {
	par(mfrow=c(1,2), mar=c(6,5,1,1))
	
	# boxplot
	boxplot(mut$dia, mut$rel, xlab="", ylab="# non-silent mutations", na.action=na.exclude, outline=F, cex.axis=1.5, xaxt="n", cex.lab=1.5, ylim=c(0,100))
	axis(1, at=c(1,2), cex.axis=1.5, labels=c("diagnosis", "relapse"))
	stripchart(list(mut$dia, mut$rel), method="jitter", vertical=T, pch=19, col=c("red", "blue"), add=T)
	print(sprintf("P-value Kruskal-Wallis: %.2g", test$p.value))
	
	# barplot
	barplot(t(bardata), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, ylab="", ylim=c(0, 100), cex.axis=1.5, cex.names=1.5, cex.lab=1.5)
	abline(h=seq(10, 100, by=10), col="gray")
	box()
	barplot(t(bardata), beside=FALSE, col=brewer.pal(6, "Set2"), las=2, ylab="# non-silent mutations", ylim=c(0, 100), cex.axis=1.5, cex.names=1.5, cex.lab=1.5, add=T)
	legend(x=1, y=95, bg="white", rev(colnames(bardata)), fill=rev(brewer.pal(3, "Set2")), cex=1.5)
}

# plot PNG
png("/mnt/projects/hdall/results/figures/supp-fig-1.png", width=2800, height=1300, pointsize=40)
make_plot()
dev.off()

# plot PDF
pdf("/mnt/projects/hdall/results/figures/supp-fig-1.pdf", width=13, height=6)
make_plot()
dev.off()
