library(GenomeGraphs)
library(biomaRt)
library(gdata)

X11.options(type="Xlib")

chr <- 12

#data("exampleData", package = "GenomeGraphs")
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

d <- read.xls("/data/christian/peter/Log2 Ratios_segments.xls")

genomeAxis <- makeGenomeAxis(add53 = F, add35 = F)
probestart <- c(1000, 1000, 1000)

chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
#chromosomes <- c("13")

pdf(file="~/peter/coverage-plots.pdf", height=4, width=10)
for (chr in chromosomes)
{
	ideog <- makeIdeogram(chromosome = chr)

	segStart <- d[d$Chromosome==paste("chr", chr, sep=""), "Start"]
	segEnd <- d[d$Chromosome==paste("chr", chr, sep=""), "End"]
	segments <- d[d$Chromosome==paste("chr", chr, sep=""), "Log2Ratio"]

	# add dummy segment to at height zero to allow plot to determine valid y-range; otherwise, nothing is plotted if there is only a single segment
	segStart <- c(0, segStart)
	segEnd <- c(0, segEnd)
	segments <- c(0, segments)
	
	seg <- makeSegmentation(segStart, segEnd, segments, dp = DisplayPars(color = "black", lwd = 2, lty = "solid"))                                  
	probe.minmax <- matrix(c(min(segments), mean(segments), max(segments)))
	cop <- makeGenericArray(intensity = probe.minmax, trackOverlay=seg, probeStart = probestart, dp = DisplayPars(size = 3, color = "white", type = "dot", ylim=c(min(segments, -2), max(segments, 2))))

	mylist <- list()
	mylist[[paste("chr", chr, sep="")]] <- ideog
	mylist[[2]] <- genomeAxis
	mylist[["Log2Ratio"]] <- cop

	gdPlot(mylist, minBase = min(segStart), maxBase = max(segEnd), labelCex = 1)
}
dev.off()