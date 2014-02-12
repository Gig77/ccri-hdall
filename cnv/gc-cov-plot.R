options(warn=1)

#library(lattice)
#library(latticeExtra)

png("all.gc-cov-plot.png", width=1800, height=10000)

par(mfrow=c(20,3))
for (p in c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "A", "B", "C", "D", "E", "X", "Y")) {
	rem <- read.delim(paste0(p, "_rem.coverage.bed"), header=F)
	rem <- rem[rem$V6 > 0,]
	dia <- read.delim(paste0(p, "_dia.coverage.bed"), header=F)
	dia <- dia[dia$V6 > 0,]
	rel <- read.delim(paste0(p, "_rel.coverage.bed"), header=F)
	rel <- rel[rel$V6 > 0,]
	
	#xyplot(rem$V6~rem$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rem"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	#xyplot(dia$V6~dia$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " dia"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	#xyplot(rel$V6~rel$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rel"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	
	plot(rem$V5, rem$V6, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rem"), log="y", ylim=c(10,10000), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
	plot(dia$V5, dia$V6, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " dia"), log="y", ylim=c(10,10000), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
	plot(rel$V5, rel$V6, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rel"), log="y", ylim=c(10,10000), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)	
}

dev.off()

