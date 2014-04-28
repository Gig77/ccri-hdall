options(warn=1)

#library(lattice)
#library(latticeExtra)

# coverage avoid bias due to aneuploidy in tumor samples 
diploid.chr <- list("430"=c("chr2", "chr3", "chr8", "chr9", "chr12", "chr13", "chr15", "chr16", "chr17", "chr19", "chr20", "chr22"))

#patients <- c("430")
#png("430.gc-cov-plot.allchr.png", width=1800, height=500)
#par(mfrow=c(1, 3))
#alpha <- 0.05
#only.disomic <- FALSE

#patients <- c("430")
#png("430.gc-cov-plot.disomic.png", width=1800, height=500)
#par(mfrow=c(1, 3))
#alpha <- 0.1
#only.disomic <- TRUE

patients <- c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "A", "B", "C", "D", "E", "X", "Y")
png("all.gc-cov-plot.png", width=1800, height=10000)
par(mfrow=c(20,3))
alpha <- 0.05
only.disomic <- FALSE

for (p in patients) {
	rem <- read.delim(paste0(p, "_rem.coverage.bed"), header=F)
	rem <- rem[rem$V6 > 0,]
	if (only.disomic && !is.null(diploid.chr[[p]])) rem <- rem[rem$V1 %in% diploid.chr[[p]],]
	dia <- read.delim(paste0(p, "_dia.coverage.bed"), header=F)
	dia <- dia[dia$V6 > 0,]
	if (only.disomic && !is.null(diploid.chr[[p]])) dia <- dia[dia$V1 %in% diploid.chr[[p]],]
	rel <- read.delim(paste0(p, "_rel.coverage.bed"), header=F)
	rel <- rel[rel$V6 > 0,]
	if (only.disomic && !is.null(diploid.chr[[p]])) rel <- rel[rel$V1 %in% diploid.chr[[p]],]
	
	#xyplot(rem$V6~rem$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rem"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	#xyplot(dia$V6~dia$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " dia"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	#xyplot(rel$V6~rel$V5, scales=list(y=list(log=10)), yscale.components=yscale.components.log10ticks, cex=0.1, col=rgb(0,0,0,0.05), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rel"), cex.lab=1.5, cex.main=1.5, ylim=c(10, 10000))
	
	plot(rem$V5, rem$V6, cex=0.1, col=rgb(0,0,0,alpha), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rem, ", ifelse(only.disomic, "disomic chromosomes", "all chromosomes")), log="y", ylim=c(10,10000), xlim=c(0, 1), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
	plot(dia$V5, dia$V6, cex=0.1, col=rgb(0,0,0,alpha), xlab="GC content", ylab="log10 coverage", main=paste0(p, " dia, ", ifelse(only.disomic, "disomic chromosomes", "all chromosomes")), log="y", ylim=c(10,10000), xlim=c(0, 1), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
	plot(rel$V5, rel$V6, cex=0.1, col=rgb(0,0,0,alpha), xlab="GC content", ylab="log10 coverage", main=paste0(p, " rel, ", ifelse(only.disomic, "disomic chromosomes", "all chromosomes")), log="y", ylim=c(10,10000), xlim=c(0, 1), cex.lab=1.5, cex.main=1.5, cex.sub=1.5)	
}

dev.off()

