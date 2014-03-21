# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'

maxdepth <- 3000

# Get a list of the bedtools output files you'd like to read in
files <- list.files(path="~/hdall/results/reseq/coverage", pattern=".coverage.bedtools.txt$")
labs <- sapply(strsplit(files, ".", fixed=T), "[[", 1) # extract sample name from file name
labs <- gsub("Diagnosis", "dia", labs)
labs <- gsub("Relapse", "rel", labs)
labs <- gsub("Remission", "rem", labs)

# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
means <- numeric(0)
for (i in 1:length(files)) {
	cov[[i]] <- read.table(paste0("~/hdall/results/reseq/coverage/", files[i]))
	cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
	means[i] <- cov_cumul[[i]][1000]
}

# Pick some colors
# Ugly:
# cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
library(RColorBrewer)
#cols <- brewer.pal(length(cov), "Dark2")
cols <- rainbow(length(cov))
ltypes <- rep(1:6,length.out=length(cov))

# Save the graph to a file
png("~/hdall/results/reseq/coverage/coverage-plots.png", h=2000, w=2500, pointsize=40)

# Create plot area, but do not plot anything. Add gridlines and axis labels.
layout(matrix(c(1,2), nrow = 1), widths = c(0.75, 0.25))
par()
plot(cov[[1]][2:(maxdepth+1), 2], cov_cumul[[1]][1:maxdepth], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Unique Coverage", xaxt="n", yaxt="n")
abline(v = c(0, 250, 500, 1000, 1500, 2000, 2500, 3000), col = "gray60", lty=3)
abline(h = seq(0, 1, 0.1), col = "gray60", lty=3)
abline(h = 0.5, col = "gray60", lty=1)
axis(1, at=c(0, 250, 500, 1000, 1500, 2000, 2500, 3000))
axis(2, at=seq(0, 1, 0.1))

# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:(maxdepth+1), 2], cov_cumul[[i]][1:maxdepth], type='l', lwd=2, lty=ltypes[i], col=cols[i])

# Add a legend using the nice sample labeles rather than the full filenames.
par(mar=c(0, 0, 8, 0), cex=0.5)
plot(0:1, 0:1, type="n", axes=F, ann=F)
legend("topleft", legend=labs[order(means, decreasing=T)], col=cols[order(means, decreasing=T)], lwd=3, lty=ltypes[order(means, decreasing=T)], ncol=2)

dev.off()