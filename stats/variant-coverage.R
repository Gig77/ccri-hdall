# TABLE: filtered-variants.cosmic.normaf.tsv
v <- read.delim("/mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv")
v <- v[v$status != "REJECT",]

dp.dia <- v[v$sample=="rem_dia", c("patient", "dp_leu_tot")]
dp.rel <- v[v$sample=="rem_rel", c("patient", "dp_leu_tot")]
dp.rem <- unique(v[,c("patient", "chr", "pos", "dp_rem_tot")])[,c("patient", "dp_rem_tot")]

pdf("/mnt/projects/hdall/results/stats/variant-coverage.pdf")
par(mfrow=c(3,1), mar=(c(2, 3, 2, 0.5)))

plot(dp.dia$patient, dp.dia$dp_leu_tot, ylim=c(0,300), yaxt="n", main=paste("dia (n=", nrow(dp.dia), ")", sep=""), cex.axis=0.8)
points(unique(as.numeric(dp.dia$patient)), tapply(dp.dia$dp_leu_tot, dp.dia$patient, mean), col="blue")
axis(2, at = seq(0, 300, 25), las = 1, cex.axis=0.8) 
abline(h=median(dp.dia$dp_leu_tot), col="red")
abline(h=mean(dp.dia$dp_leu_tot), col="blue")

plot(dp.rel$patient, dp.rel$dp_leu_tot, ylim=c(0,300), yaxt="n", main=paste("rel (n=", nrow(dp.rel), ")", sep=""), cex.axis=0.8)
points(unique(as.numeric(dp.rel$patient)), tapply(dp.rel$dp_leu_tot, dp.rel$patient, mean), col="blue")
axis(2, at = seq(0, 300, 25), las = 1, cex.axis=0.8) 
abline(h=median(dp.rel$dp_leu_tot), col="red")
abline(h=mean(dp.rel$dp_leu_tot), col="blue")

plot(dp.rem$patient, dp.rem$dp_rem_tot, ylim=c(0,300), yaxt="n", main=paste("rem (n=", nrow(dp.rem), ")", sep=""), cex.axis=0.8)
points(unique(as.numeric(dp.rem$patient)), tapply(dp.rem$dp_rem_tot, dp.rem$patient, mean), col="blue")
axis(2, at = seq(0, 300, 25), las = 1, cex.axis=0.8) 
abline(h=median(dp.rem$dp_rem_tot), col="red")
abline(h=mean(dp.rem$dp_rem_tot), col="blue")

dev.off()