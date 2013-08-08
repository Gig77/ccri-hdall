options(warn=1)

chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
#sv <- read.delim("~/hdall/results/filtered-variants.cosmic.normaf.tsv")

pdf("~/hdall/results/stats/variants-per-chrom-and-ploidy.pdf", height=10)

for (patient in c("314", "399", "430", "446", "460", "545", "592", "715", "818", "B", "C"))
{
	par(mfcol=c(2,1))
	for (sample in c("dia", "rel"))
	{
		s <- paste("rem_", sample, sep="")
		sv.chr.monosom <- table(sv[sv$sample==s & sv$copy_no==1,c("patient", "chr")])
		sv.chr.disom <- table(sv[sv$sample==s & sv$copy_no==2,c("patient", "chr")])
		sv.chr.trisom <- table(sv[sv$sample==s & sv$copy_no==3,c("patient", "chr")])
		sv.chr.tetrasom <- table(sv[sv$sample==s & sv$copy_no==4,c("patient", "chr")])
		n <- sum(sv$patient==patient & sv$sample==s)
		
		barplot(
			rbind(
				sv.chr.monosom[patient,chr], 
				sv.chr.disom[patient,chr], 
				sv.chr.trisom[patient,chr], 
				sv.chr.tetrasom[patient,chr]
			), 
			col=c("black", "gray", "red", "blue"),
			las=2,
			ylab="# somatic variants", 
			main=paste(patient, " ", sample, " (n=", n, ")", sep=""), 
			legend.text=c(
				paste("monosom (", sum(sv$patient==patient & sv$sample==s & sv$copy_no==1), ")", sep=""), 
				paste("disom (", sum(sv$patient==patient & sv$sample==s & sv$copy_no==2), ")", sep=""), 
				paste("trisom (", sum(sv$patient==patient & sv$sample==s & sv$copy_no==3), ")", sep=""), 
				paste("tetrasom (", sum(sv$patient==patient & sv$sample==s & sv$copy_no==4), ")", sep="")
			) 
		)
	}
}

dev.off()