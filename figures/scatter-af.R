options(warn=1)

mincof = 10

pdf("~/hdall/results/reseq/scatter-af.pdf")
#png("~/hdall/results/reseq/scatter-af.png", h=2000, w=2000, pointsize=40)
m <- read.delim("~/hdall/results/reseq/filtered-variants.merged.tsv", check.names=F, stringsAsFactor=F)
mf <- m[!is.na(m$status.dis) & m$status.dis=="PASS" & !is.na(m$status.val) & m$status.val=="PASS",]
mf <- mf[mf$dp_leu_tot.dis >= mincof & mf$dp_leu_tot.val >= mincof,]
fit <-lm(mf$freq_leu.dis~mf$freq_leu.val)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(mf$freq_leu.dis, mf$freq_leu.val, ylim=c(0, 1), xlim=c(0, 1), ylab="AF exome sequencing", xlab="AF targeted sequencing", main=sprintf("R=%.2f, p=%.2g", R, p))
abline(fit, col="red")
dev.off()
