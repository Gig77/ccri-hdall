options(warn=1)

library(MASS)

# patient 545
patient <- "545";
sample <- "dia";
diploid <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr7", "chr8", "chr9", "chr11", "chr12", "chr13", "chr15", "chr16")
triploid <- c("chr6", "chr10")

#v <- read.csv(paste("~/hdall/data/mutect_vcf/", patient, "_rem_", sample, "_call_stats.out", sep=""), sep="\t", skip=1)
vpass <- v[v$t_ref_count+v$t_alt_count >= 50 & v$n_ref_count+v$n_alt_count >= 50 & v$t_ins_count == 0 & v$t_del_count == 0 & v$t_ref_max_mapq >= 60 & v$t_alt_max_mapq >= 60 , c("contig", "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
vgerm <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count > 1,]
vsom <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count == 0,]

# triploid chromosomes
#----------------------
vtrip <- vgerm[vgerm$contig %in% triploid,]
af <- vtrip$t_alt_count / (vtrip$t_ref_count + vtrip$t_alt_count)
afhet.trip <- af[af>=0.23 & af<=0.43]	

# fit normal distribution
f.trip<-fitdistr(afhet.trip, "normal")
x <- seq(0,1,length=300)
y <- dnorm(x, f.trip$estimate["mean"]-0.05, f.trip$estimate["sd"])
hist(af, xlim=c(0, 1), xlab="allelic frequency", breaks=50, main=paste("patient ", patient, " ", sample))
axis(1, at=seq(0,1,0.1)) 
par(new=T)
plot(x, y, type="l", lwd=2, col="red", axes=F, xlab=NA, ylab=NA)

freqtable <- cbind(rep(patient, 100), rep(sample, 100), rep("triploid", 100), seq(0.01,1,0.01), pnorm(seq(0.01,1,0.01), f.trip$estimate["mean"], f.trip$estimate["sd"]))
colnames(freqtable) <- c("patient", "sample", "ploidy", "af", "probability")
write.table(freqtable, file=paste("~/hdall/results/clonal-analysis/allelic-freq-prob.tsv", sep=""), col.names=T, row.names=F, sep="\t", quote=F)

stop()

# diploid chromosomes
#----------------------
vdipl <- vgerm[vgerm$contig %in% diploid,]
af <- vdipl$t_alt_count / (vdipl$t_ref_count + vdipl$t_alt_count)
afhet.dip <- af[af>=0.3 & af<=0.7]	

# determine most frequent allelic frequency assuming that this corresponds to the peak of heterozygous variants
peak <- as.numeric( sub("\\((.+),.*", "\\1", names(table(cut(af, 100)))[50]))

# get allelic frequency above peak and below homozygous variants
afhet.dip.onetailed <- af[af>=peak & af<=0.8]	

af.freq.onetailed <- table(cut(afhet.dip.onetailed, breaks=30))
stepsize <- as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", names(af.freq.onetailed)[1])) - as.numeric(sub("\\((.+),.*", "\\1", names(af.freq.onetailed)[1]))

# symmetrify distribution and shift bin to left for optimal fit
af.freq.twotailed <- c(rep(0,19),rev(af.freq.onetailed), af.freq.onetailed,rep(0,21))

af.cumfreq <- cumsum(af.freq.twotailed)/length(afhet.dip.onetailed)/2
freqtable <- cbind(seq(0.01,1,0.01),af.cumfreq)
colnames(freqtable) <- c("af", "probability")
write.table(freqtable, file=paste("~/hdall/results/clonal-analysis/allelic-freq-prob.", patient, ".disomic.", sample, ".tsv", sep=""), col.names=T, row.names=F, sep="\t", quote=F)
#stop()

pdf(paste("~/hdall/results/clonal-analysis/allelic-freq-prob.", patient, ".disomic.", sample, ".distribution.pdf", sep=""))

# plot histogram
hist(af, xlim=c(0, 1), xlab="allelic frequency", breaks=50, main=paste("patient ", patient, " ", sample))
axis(1, at=seq(0,1,0.1)) 

# fit normal distribution
f<-fitdistr(afhet.dip, "normal")
x <- seq(0,1,length=300)
y <- dnorm(x, f$estimate["mean"], f$estimate["sd"])
par(new=T)
plot(x, y, type="l", lwd=2, col="red", axes=F, xlab=NA, ylab=NA)

# fit fat-tailed cauchy distribution
f<-fitdistr(afhet.dip, "cauchy")
x <- seq(0,1,length=300)
y <- dcauchy(x, f$estimate["location"], f$estimate["scale"])
par(new=T)
plot(x, y, type="l", lwd=2, col="blue", axes=F, xlab=NA, ylab=NA)

#par(new=T)
#plot(seq(0.01,1,0.01), af.freq.twotailed, axes=F, xlab=NA, ylab=NA)

dev.off()