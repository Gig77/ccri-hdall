# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'

maxdepth.exome <- 200
maxdepth.panel <- 1500

# Get a list of the bedtools output files you'd like to read in
files.exome <- list.files(path="/mnt/projects/hdall/results/coverage", pattern=".coverage.bedtools.txt$")
files.exome <- files.exome[! files.exome %in% c("E_rem.coverage.bedtools.txt", "E_dia.coverage.bedtools.txt", "E_rel.coverage.bedtools.txt")] # remove patient E
labs.exome <- sapply(strsplit(files.exome, ".", fixed=T), "[[", 1) # extract sample name from file name

files.panel.relapsing <- list.files(path="/mnt/projects/hdall/results/reseq/coverage", pattern=".coverage.bedtools.txt$")
labs.panel <- sapply(strsplit(files.panel.relapsing, ".", fixed=T), "[[", 1) # extract sample name from file name
files.panel.nonrelapsing <- list.files(path="/mnt/projects/hdall/results/reseq/coverage/nonrelapsing", pattern=".coverage.bedtools.txt$")
labs.panel <- c(labs.panel, sapply(strsplit(files.panel.relapsing, ".", fixed=T), "[[", 1)) # extract sample name from file name
labs.panel <- gsub("Diagnosis", "dia", labs.panel)
labs.panel <- gsub("Relapse", "rel", labs.panel)
labs.panel <- gsub("Remission", "rem", labs.panel)

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov.exome <- list()
cov_cumul.exome <- list()
means.exome <- numeric(0)
for (i in 1:length(files.exome)) {
	cov.exome[[i]] <- read.table(paste0("/mnt/projects/hdall/results/coverage/", files.exome[i]))
	cov_cumul.exome[[i]] <- 1-cumsum(cov.exome[[i]][,5])
	means.exome[i] <- cov_cumul.exome[[i]][50]
}

cov.panel <- list()
cov_cumul.panel <- list()
means.panel <- numeric(0)
for (i in 1:length(files.panel.relapsing)) {
	cov.panel[[i]] <- read.table(paste0("/mnt/projects/hdall/results/reseq/coverage/", files.panel.relapsing[i]))
	cov_cumul.panel[[i]] <- 1-cumsum(cov.panel[[i]][,5])
	means.panel[i] <- cov_cumul.panel[[i]][250]
}
for (i in 1:length(files.panel.nonrelapsing)) {
	idx <- i + length(files.panel.relapsing)
	cov.panel[[idx]] <- read.table(paste0("/mnt/projects/hdall/results/reseq/coverage/nonrelapsing/", files.panel.nonrelapsing[i]))
	cov_cumul.panel[[idx]] <- 1-cumsum(cov.panel[[idx]][,5])
	means.panel[idx] <- cov_cumul.panel[[idx]][250]
}

# Pick some colors
# Ugly:
# cols <- 1:length(cov.exome)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
library(RColorBrewer)
#cols <- brewer.pal(length(cov.exome), "Dark2")
cols.exome <- rainbow(length(cov.exome))
ltypes.exome <- rep(1:6,length.out=length(cov.exome))
cols.panel <- rainbow(length(cov.panel))
ltypes.panel <- rep(1:6,length.out=length(cov.panel))

do_plot <- function(lwidth) {
	# Create plot area, but do not plot anything. Add gridlines and axis labels.
	layout(matrix(c(1,2,3,4), nrow = 1), widths = c(0.39, 0.07, 0.39, 0.15))
	par(mar=c(5.1, 4.1, 4.1, 2.1), cex=1.0)
	plot(cov.exome[[1]][2:(maxdepth.exome+1), 2], cov_cumul.exome[[1]][1:maxdepth.exome], type='n', xlab="Unique coverage", ylab="Percentage of target bases >= coverage", ylim=c(0,1.0), main="Exome sequencing", xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5)
	abline(v = c(0, 10, 25, 50, 100, 200), col = "gray60", lty=3)
	abline(h = seq(0, 1, 0.1), col = "gray60", lty=3)
	abline(h = 0.5, col = "gray60", lty=1)
	axis(1, at=c(0, 10, 25, 50, 100, 200), cex.axis=1.3)
	axis(2, at=seq(0, 1, 0.1), cex.axis=1.3)
	
	# Actually plot the data for each of the alignments (stored in the lists).
	for (i in 1:length(cov.exome)) points(cov.exome[[i]][2:(maxdepth.exome+1), 2], cov_cumul.exome[[i]][1:maxdepth.exome], type='l', lwd=lwidth, lty=ltypes.exome[i], col=cols.exome[i])
	
	# Add a legend using the nice sample labeles rather than the full filenames.
	par(mar=c(0, 0, 6, 0), cex=0.7)
	plot(0:1, 0:1, type="n", axes=F, ann=F)
	legend("topleft", legend=labs.exome[order(means.exome, decreasing=T)], col=cols.exome[order(means.exome, decreasing=T)], lwd=lwidth, lty=ltypes.exome[order(means.exome, decreasing=T)], ncol=1)	

	par(mar=c(5.1, 4.1, 4.1, 2.1), cex=1.0)
	plot(cov.panel[[1]][2:(maxdepth.panel+1), 2], cov_cumul.panel[[1]][1:maxdepth.panel], type='n', xlab="Unique coverage", ylab="Percentage of target bases >= coverage", ylim=c(0,1.0), main="Gene panel", xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5)
	abline(v = c(0, 100, 250, 500, 1000, 1500), col = "gray60", lty=3)
	abline(h = seq(0, 1, 0.1), col = "gray60", lty=3)
	abline(h = 0.5, col = "gray60", lty=1)
	axis(1, at=c(0, 100, 250, 500, 1000, 1500), cex.axis=1.3)
	axis(2, at=seq(0, 1, 0.1), cex.axis=1.3)
	
	# Actually plot the data for each of the alignments (stored in the lists).
	for (i in 1:length(cov.panel)) points(cov.panel[[i]][2:(maxdepth.panel+1), 2], cov_cumul.panel[[i]][1:maxdepth.panel], type='l', lwd=lwidth, lty=ltypes.panel[i], col=cols.panel[i])
	
	# Add a legend using the nice sample labeles rather than the full filenames.
	par(mar=c(0, 0, 8, 0), cex=0.5)
	plot(0:1, 0:1, type="n", axes=F, ann=F)
	legend("topleft", legend=labs.panel[order(means.panel, decreasing=T)], col=cols.panel[order(means.panel, decreasing=T)], lwd=lwidth, lty=ltypes.panel[order(means.panel, decreasing=T)], ncol=4)
	
}

# Save the graph to a file
png("/mnt/projects/hdall/results/coverage/coverage-plots-exome-and-panel.png", h=2500, w=6000, pointsize=40)
do_plot(3)
dev.off()

pdf("/mnt/projects/hdall/results/coverage/coverage-plots-exome-and-panel.pdf", width=6.3, height=2.6)

	# Create plot area, but do not plot anything. Add gridlines and axis labels.
	layout(matrix(c(1,2,3,4), nrow = 1), widths = c(0.38, 0.07, 0.38, 0.17))
	par(mar=c(5.1, 4.1, 4.1, 2.1), cex=0.3)
	plot(cov.exome[[1]][2:(maxdepth.exome+1), 2], cov_cumul.exome[[1]][1:maxdepth.exome], type='n', xlab="Unique coverage", ylab="Percentage of target bases >= coverage", ylim=c(0,1.0), main="Exome sequencing", xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5, bty="n")
	box(lwd=0.5)
	abline(v = c(0, 10, 25, 50, 100, 200), col = "gray60", lty=3, lwd=0.5)
	abline(h = seq(0, 1, 0.1), col = "gray60", lty=3, lwd=0.5)
	abline(h = 0.5, col = "gray60", lty=1, lwd=0.5)
	axis(1, at=c(0, 10, 25, 50, 100, 200), cex.axis=1.3, lwd=0.5)
	axis(2, at=seq(0, 1, 0.1), cex.axis=1.3, lwd=0.5)
	
	# Actually plot the data for each of the alignments (stored in the lists).
	for (i in 1:length(cov.exome)) points(cov.exome[[i]][2:(maxdepth.exome+1), 2], cov_cumul.exome[[i]][1:maxdepth.exome], type='l', lwd=0.6, lty=ltypes.exome[i], col=cols.exome[i])
	
	# Add a legend using the nice sample labeles rather than the full filenames.
	par(mar=c(0, 0, 9.3, 0), cex=0.13)
	plot(0:1, 0:1, type="n", axes=F, ann=F, bty="n")
	legend("topleft", legend=labs.exome[order(means.exome, decreasing=T)], col=cols.exome[order(means.exome, decreasing=T)], lwd=0.6, lty=ltypes.exome[order(means.exome, decreasing=T)], ncol=1, box.lwd=0.5)	
	
	par(mar=c(5.1, 4.1, 4.1, 2.1), cex=0.3)
	plot(cov.panel[[1]][2:(maxdepth.panel+1), 2], cov_cumul.panel[[1]][1:maxdepth.panel], type='n', xlab="Unique coverage", ylab="Percentage of target bases >= coverage", ylim=c(0,1.0), main="Gene panel", xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.5, bty="n")
	box(lwd=0.5)
	abline(v = c(0, 100, 250, 500, 1000, 1500), col = "gray60", lty=3, lwd=0.5)
	abline(h = seq(0, 1, 0.1), col = "gray60", lty=3, lwd=0.5)
	abline(h = 0.5, col = "gray60", lty=1, lwd=0.5)
	axis(1, at=c(0, 100, 250, 500, 1000, 1500), cex.axis=1.3, lwd=0.5)
	axis(2, at=seq(0, 1, 0.1), cex.axis=1.3, lwd=0.5)
	
	# Actually plot the data for each of the alignments (stored in the lists).
	for (i in 1:length(cov.panel)) points(cov.panel[[i]][2:(maxdepth.panel+1), 2], cov_cumul.panel[[i]][1:maxdepth.panel], type='l', lwd=0.6, lty=ltypes.panel[i], col=cols.panel[i])
	
	# Add a legend using the nice sample labeles rather than the full filenames.
	par(mar=c(0, 0, 9.3, 0), cex=0.13)
	plot(0:1, 0:1, type="n", axes=F, ann=F, bty="n")
	legend("topleft", legend=labs.panel[order(means.panel, decreasing=T)], col=cols.panel[order(means.panel, decreasing=T)], lwd=0.6, lty=ltypes.panel[order(means.panel, decreasing=T)], ncol=4, box.lwd=0.5)

dev.off()
