options(warn=1)

library(MASS)

ploidy <- data.frame(patient=character(), ploidy=character(), chr=character(), stringsAsFactors=F)

# patient 545
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr1")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr2")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr3")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr4")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr5")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr7")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr8")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr9")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr11")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr12")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr13")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr15")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr16")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr17")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr19")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr20")
ploidy[nrow(ploidy)+1,] <- c("545", "diploid", "chr22")
ploidy[nrow(ploidy)+1,] <- c("545", "triploid", "chr6")
ploidy[nrow(ploidy)+1,] <- c("545", "triploid", "chr10")

# patient 430
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr2")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr3")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr8")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr9")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr12")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr13")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr15")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr16")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr17")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr19")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr20")
ploidy[nrow(ploidy)+1,] <- c("430", "diploid", "chr22")
ploidy[nrow(ploidy)+1,] <- c("430", "triploid", "chr4")
ploidy[nrow(ploidy)+1,] <- c("430", "triploid", "chr6")
ploidy[nrow(ploidy)+1,] <- c("430", "triploid", "chr10")

# patient 446
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr2")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr3")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr13")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr15")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr19")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr20")
ploidy[nrow(ploidy)+1,] <- c("446", "diploid", "chr22")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr4")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr6")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr8")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr11")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr12")
ploidy[nrow(ploidy)+1,] <- c("446", "triploid", "chr14")

pdf("/mnt/projects/hdall/results/clonal-analysis/allelic-freq-prob.distributions.pdf")

iteration <- 1
for (patient in unique(ploidy$patient))
{
	for (sample in c("dia", "rel")) 
	{
		v <- read.csv(paste("/mnt/projects/hdall/data/mutect_vcf/", patient, "_rem_", sample, "_call_stats.out", sep=""), sep="\t", skip=1)
		vpass <- v[v$t_ref_count+v$t_alt_count >= 30 & v$n_ref_count+v$n_alt_count >= 30 & v$t_ins_count == 0 & v$t_del_count == 0 & v$t_ref_max_mapq >= 60 & v$t_alt_max_mapq >= 60 , c("contig", "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
		vgerm <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count > 1,]
		vsom <- vpass[vpass$t_alt_count > 1 & vpass$n_alt_count == 0,]
		
		
		# diploid chromosomes
		#----------------------
		diploid <- ploidy[ploidy$patient==patient & ploidy$ploidy=="diploid","chr"]
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
		write.table(freqtable, file=paste("/mnt/projects/hdall/results/clonal-analysis/allelic-freq-prob.tsv", sep=""), col.names=(iteration==1), row.names=F, sep="\t", quote=F, append=(iteration>1))
		
		# triploid chromosomes
		#----------------------
		triploid <- ploidy[ploidy$patient==patient & ploidy$ploidy=="triploid","chr"]
		vtrip.germ <- vgerm[vgerm$contig %in% triploid,]
		vtrip.som <- vsom[vsom$contig %in% triploid,]
		af.germ.trip <- vtrip.germ$t_alt_count / (vtrip.germ$t_ref_count + vtrip.germ$t_alt_count)
		af.germ.trip.pruned <- af.germ.trip[af.germ.trip>=0.20 & af.germ.trip<=0.46]	
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
		write.table(freqtable, file=paste("/mnt/projects/hdall/results/clonal-analysis/allelic-freq-prob.tsv", sep=""), col.names=F, row.names=F, sep="\t", quote=F, append=T)
		
		iteration <- iteration + 1
	}
}

dev.off()