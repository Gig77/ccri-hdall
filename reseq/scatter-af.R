options(warn=1)

pdf("scatter-af.pdf")
m <- read.delim("filtered-variants.merged.tsv", check.names=F, stringsAsFactor=F)
mf <- m[!is.na(m$status.dis) & m$status.dis=="PASS" & !is.na(m$status.val) & m$status.val=="PASS",]
plot(mf$freq_leu.dis, mf$freq_leu.val, ylim=c(0, 1), xlim=c(0, 1), ylab="AF exome sequencing", xlab="AF targeted sequencing", main=paste("R=", cor(mf$freq_leu.dis, mf$freq_leu.val), sep=""))
reg <-lm(mf$freq_leu.dis~mf$freq_leu.val)
abline(reg, col="red")
dev.off()
