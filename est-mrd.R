options(warn=1)

# estimate minimal residual disease
patients <- c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "1021247", "A", "B", "C", "D", "E", "X", "Y") 
#patients <- c("399", "314")

pdf("~/hdall/results/mrd-rel.pdf", width=12, paper='A4r')
par(mfrow=c(4,5), mar=c(2,1,2,1))
for(p in patients) {
	t <- read.csv(paste("~/hdall/data/mutect_vcf/", p, "_rem_rel_call_stats.out", sep=""), sep="\t", skip=1)
	tmajor <- t[t$t_ref_count+t$t_alt_count >= 50 & t$n_ref_count+t$n_alt_count >= 50 & t$t_ins_count == 0 & t$t_del_count == 0 & t$t_ref_max_mapq >= 60 & t$t_alt_max_mapq >= 60 , c("t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
	frac_norm <- tmajor$n_alt_count/(tmajor$n_ref_count+tmajor$n_alt_count)
	frac_tumor <- tmajor$t_alt_count/(tmajor$t_ref_count+tmajor$t_alt_count)
	smoothScatter(frac_norm, frac_tumor, xlim=c(0,1), ylim=c(0,1), main=paste(p, " (n=", length(frac_norm), ")", sep=""))
		
#	nfrac <- nfrac[!is.na(nfrac) & nfrac < 0.25]
#	nfrac <- data.frame(rep(p, length(nfrac)), nfrac) 
#	plot(density(nfrac), main=paste(p, " (n=", length(nfrac), ")", sep=""))	
}
dev.off()

pdf("~/hdall/results/mrd-dia.pdf", width=12, paper='A4r')
par(mfrow=c(4,5), mar=c(2,0.5,2,0.5))
for(p in patients) {
	t <- read.csv(paste("~/hdall/data/mutect_vcf/", p, "_rem_dia_call_stats.out", sep=""), sep="\t", skip=1)
	tmajor <- t[t$t_ref_count+t$t_alt_count >= 50 & t$n_ref_count+t$n_alt_count >= 50 & t$t_ins_count == 0 & t$t_del_count == 0 & t$t_ref_max_mapq >= 60 & t$t_alt_max_mapq >= 60 , c("t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
	frac_norm <- tmajor$n_alt_count/(tmajor$n_ref_count+tmajor$n_alt_count)
	frac_tumor <- tmajor$t_alt_count/(tmajor$t_ref_count+tmajor$t_alt_count)
	smoothScatter(frac_norm, frac_tumor, xlim=c(0,1), ylim=c(0,1), main=paste(p, " (n=", length(frac_norm), ")", sep=""))
		
#	nfrac <- nfrac[!is.na(nfrac) & nfrac < 0.25]
#	nfrac <- data.frame(rep(p, length(nfrac)), nfrac) 
#	plot(density(nfrac), main=paste(p, " (n=", length(nfrac), ")", sep=""))	
}
dev.off()
