options(warn=1)

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("patient not specified")
if (is.na(args[2])) stop("diagnosis sample not specified")
if (is.na(args[3])) stop("relapse sample not specified")
if (is.na(args[4])) stop("remission sample not specified")
patient <- args[1]

dia <- read.delim(args[2], header=F)
rel <- read.delim(args[3], header=F)
rem <- read.delim(args[4], header=F)

names(dia) <- c("chr", "start", "end", "name", "length", "strand", "total.dia", "avg.dia")
names(rem) <- c("chr", "start", "end", "name", "length", "strand", "total.rem", "avg.rem")
names(rel) <- c("chr", "start", "end", "name", "length", "strand", "total.rel", "avg.rel")

m <- merge(rem, dia, by=c("chr", "start", "end", "name", "length", "strand"))
m <- merge(m, rel, by=c("chr", "start", "end", "name", "length", "strand"))
m <- cbind(m, sapply(strsplit(as.character(m$name), "-"), "[", 1))
m <- cbind(m, log(m$avg.dia/m$avg.rem, 2))
m <- cbind(m, log(m$avg.rel/m$avg.rem, 2))
m <- m[with(m, order(chr, start)),]

names(m)[13] <- "gene"
names(m)[14] <- "log.rem.dia"
names(m)[15] <- "log.rem.rel"

pdf(paste("cov-plot.", patient, ".reseq.pdf.part", sep=""), width=15)
par(mfrow=c(2,1))
qplot(gene, log.rem.dia, data=m) + ylim(-1.5, 1.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
qplot(gene, log.rem.rel, data=m) + ylim(-1.5, 1.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

### below code is copy-paste from R console; needs to be cross-checked against the above code and cleaned up

head(m)
cbind(m, log(m$avg.dia/m$avg.rem))
cbind(m, log(m$avg.dia/m$avg.rem, 2))
m <- cbind(m, log(m$avg.dia/m$avg.rem, 2))
head(m)
as.factor(sapply(strsplit(as.character(m$name), "-"), "[", 1))
cbind(m, as.factor(sapply(strsplit(as.character(m$name), "-"), "[", 1)))
m <- cbind(m, as.factor(sapply(strsplit(as.character(m$name), "-"), "[", 1)))
head(m)
names(m)[1]
names(m)[14] <- "log.dia.rem"
head(m)
names(m)[13] <- "log.dia.rem"
names(m)[14] <- "gene"
head(m)
plot(m$gene, m$log.dia.rem)
X11.options(type="Xlib")
plot(m$gene, m$log.dia.rem)
xyplot(m$gene, m$log.dia.rem)
require(ggplot2)
install.packages("ggplot2")
library(ggplot2)
qplot(m$gene, m$log.dia.rem)
?qplot
qplot(m$gene, m$log.dia.rem) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
history()
