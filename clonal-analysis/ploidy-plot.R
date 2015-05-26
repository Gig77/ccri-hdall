options(warn=1)

# estimate minimal residual disease
patients <- c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "1021247", "A", "B", "C", "D", "E", "X", "Y") 
#patients <- c("430")

# per chromosome
pdf("/mnt/projects/hdall/results/clonal-analysis/ploidy-lowres.pdf", width=12, paper='A4r')
for(p in patients) {
	for(s in c("dia", "rel")) {
		par(mfrow=c(5,5), mar=c(2,2,1.5,0.5), oma=c(2.5,2.5,2.5,0))
		t <- read.csv(paste("/mnt/projects/hdall/data/mutect_vcf/", p, "_rem_", s, "_call_stats.out", sep=""), sep="\t", skip=1)
		tpass <- t[t$t_ref_count+t$t_alt_count >= 50 & t$n_ref_count+t$n_alt_count >= 50 & t$t_ins_count == 0 & t$t_del_count == 0 & t$t_ref_max_mapq >= 60 & t$t_alt_max_mapq >= 60 , c("contig", "t_ref_count", "t_alt_count", "n_ref_count", "n_alt_count")]
		tot_norm = sum(tpass$n_ref_count+tpass$n_alt_count) 
		tot_tum = sum(tpass$t_ref_count+tpass$t_alt_count)
	
		for (chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {
			tmajor <- tpass[tpass$contig == chr,]
			dp <- log2(((tmajor$t_alt_count+tmajor$t_ref_count)/tot_tum)/((tmajor$n_ref_count+tmajor$n_alt_count)/tot_norm))
			frac_tumor <- tmajor$t_alt_count/(tmajor$t_ref_count+tmajor$t_alt_count)	
			title <- paste(chr, " (n=", length(frac_tumor), ")", sep="")
			if (length(dp) > 1) { 
				smoothScatter(dp, frac_tumor, xlim=c(-2,2), nbin=50, ylim=c(0,1), axes=F, main=title)
			}
			else {
				plot(1, type="n", xlim=c(-2,2), ylim=c(0,1), axes=F, main=title)
				box()
			}
			axis(2, at=seq(0,1,0.1)) 
			axis(1, at=seq(-2,2,0.5)) 
		}
		mtext(paste(p, s), side=3, outer=T, cex=1, line=1)
		mtext("log2(coverage tumor/coverage remission)", side=1, outer=T, cex=1, line=1)
		mtext("tumor allelic frequency", side=2, outer=T, cex=1, line=1)
		frame()
	}
}
dev.off()
