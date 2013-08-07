options(warn=1)

library(MASS)

# patient 545
patient <- "545";
diploid <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr7", "chr8", "chr9", "chr11", "chr12", "chr13", "chr15", "chr16", "chr17", "chr19", "chr20", "chr22")
triploid <- c("chr6", "chr10")

pdf("~/hdall/results/clonal-analysis/allelic-freq-prob.distributions.pdf")

iteration <- 1
for (sample in c("dia", "rel")) 
{
	v <- read.csv(paste("~/hdall/data/mutect_vcf/", patient, "_rem_", sample, "_call_stats.out", sep=""), sep="\t", skip=1)
	vpass <- v[v$t_ref_count+v$t_alt_count >= 30 & v$n_ref_count+v$n_alt_count >= 30 & v$t_ins_count == 0 & v$t_del_count == 0 & v$t_ref_max_mapq >= 60 & v$t_alt_max_mapq >= 60 , c("contig", "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
	vgerm <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count > 1,]
	vsom <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count == 0,]
	
	
	# diploid chromosomes
	#----------------------
	vdipl.germ <- vgerm[vgerm$contig %in% diploid,]
	vdipl.som <- vsom[vsom$contig %in% diploid,]
	af.germ.dip <- vdipl.germ$t_alt_count / (vdipl.germ$t_ref_count + vdipl.germ$t_alt_count)
	af.som.dip <- vdipl.som$t_alt_count / (vdipl.som$t_ref_count + vdipl.som$t_alt_count)
	af.germ.dip.pruned <- af.germ.dip[af.germ.dip>=0.3 & af.germ.dip<=0.7]	
	af.som.dip.pruned <- af.som.dip[af.som.dip>=0.3 & af.som.dip<=0.7]	
	
	# plot histograms
	hist(af.germ.dip, xlim=c(0, 1), xlab="allelic frequency", breaks=seq(0,1,0.02), col=rgb(1,0,0,1/4), main=paste("patient", patient, sample, "diploid"), freq=T)
	par(new=T)
	hist(af.som.dip, xlim=c(0, 1), breaks=seq(0,1,0.02), axes=F, xlab=NA, ylab=NA, col=rgb(0,0,1,1/4), main=NA, freq=T)
	axis(1, at=seq(0,1,0.1)) 
	axis(4) 
	
	# fit and plot normal distribution for germline mutations
	f.germ.dip<-fitdistr(af.germ.dip.pruned, "normal")
	par(new=T)
	plot(seq(0,1,length=300), dnorm(seq(0,1,length=300), f.germ.dip$estimate["mean"], f.germ.dip$estimate["sd"]), type="l", lwd=2, col="red", axes=F, xlab=NA, ylab=NA)
	
	# fit and plot normal distribution for germline mutations
	f.som.dip<-fitdistr(af.som.dip.pruned, "normal")
	par(new=T)
	plot(seq(0,1,length=300), dnorm(seq(0,1,length=300), f.som.dip$estimate["mean"], f.germ.dip$estimate["sd"]), type="l", lwd=2, col="blue", axes=F, xlab=NA, ylab=NA)
	
	# add legend
	legend("topright", c("germline", "somatic"), lty=1, col = c("red", "blue"), merge = TRUE, inset=0.05)
	
	# output probabilities for mean-shifted germline distribution
	freqtable <- cbind(rep(patient, 100), rep(sample, 100), rep("diploid", 100), seq(0.01,1,0.01), pnorm(seq(0.01,1,0.01), f.som.dip$estimate["mean"], f.germ.dip$estimate["sd"]))
	colnames(freqtable) <- c("patient", "sample", "ploidy", "af", "probability")
	write.table(freqtable, file=paste("~/hdall/results/clonal-analysis/allelic-freq-prob.tsv", sep=""), col.names=(iteration==1), row.names=F, sep="\t", quote=F, append=(iteration>1))
	
	# triploid chromosomes
	#----------------------
	vtrip.germ <- vgerm[vgerm$contig %in% triploid,]
	vtrip.som <- vsom[vsom$contig %in% triploid,]
	af.germ.trip <- vtrip.germ$t_alt_count / (vtrip.germ$t_ref_count + vtrip.germ$t_alt_count)
	af.germ.trip.pruned <- af[af>=0.20 & af<=0.46]	
	af.som.trip <- vtrip.som$t_alt_count / (vtrip.som$t_ref_count + vtrip.som$t_alt_count)
	
	# plot histograms
	hist(af.germ.trip, xlim=c(0, 1), xlab="allelic frequency", breaks=seq(0,1,0.02), col=rgb(1,0,0,1/4), main=paste("patient", patient, sample, "triploid"))
	par(new=T)
	hist(af.som.trip, xlim=c(0, 1), breaks=seq(0,1,0.02), axes=F, xlab=NA, ylab=NA, col=rgb(0,0,1,1/4), main=NA)
	axis(1, at=seq(0,1,0.1)) 
	axis(4) 
	
	# fit and plot normal distribution for germline mutations
	f.germ.trip <- fitdistr(af.germ.trip.pruned, "normal")
	par(new=T)
	plot(seq(0,1,length=300), dnorm(seq(0,1,length=300), f.germ.trip$estimate["mean"], f.germ.trip$estimate["sd"]), type="l", lwd=2, col="red", axes=F, xlab=NA, ylab=NA)
	par(new=T)
	meanshift <- f.germ.dip$estimate["mean"] - f.som.dip$estimate["mean"]
	plot(seq(0,1,length=300), dnorm(seq(0,1,length=300), f.germ.trip$estimate["mean"]-meanshift, f.germ.trip$estimate["sd"]), type="l", lwd=2, col="blue", axes=F, xlab=NA, ylab=NA)
	
	# add legend
	legend("topright", c("germline", "somatic"), lty=1, col = c("red", "blue"), merge = TRUE, inset=0.05)
	
	freqtable <- cbind(rep(patient, 100), rep(sample, 100), rep("triploid", 100), seq(0.01,1,0.01), pnorm(seq(0.01,1,0.01), f.germ.trip$estimate["mean"]-meanshift, f.germ.trip$estimate["sd"]))
	colnames(freqtable) <- c("patient", "sample", "ploidy", "af", "probability")
	write.table(freqtable, file=paste("~/hdall/results/clonal-analysis/allelic-freq-prob.tsv", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
	
	iteration <- iteration + 1
}

dev.off()