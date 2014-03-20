options(warn=1)
library(DNAcopy)
library(optparse)

option_list <- list(
		make_option("--patient", type="character", help="patient ID"),
		make_option("--tumor", type="character", help="tumor coverage data file"),
		make_option("--normal", type="character", help="remission coverage data file"),
		make_option("--output", type="character", help="output file name"),
		make_option("--circos", type="character", help="circos output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.na(opt$patient)) stop("patient ID not specified")
if (is.na(opt$tumor)) stop("tumor coverage file not specified")
if (is.na(opt$normal)) stop("normal coverage file not specified")
if (is.na(opt$output)) stop("output file not specified")
if (is.na(opt$circos)) stop("circos output file not specified")

# for test purposes
#opt <- data.frame(patient="1020540_dia", tumor="../reseq/cnv/segmented_coverage/1020540_Diagnosis.segmented-coverage.tsv", normal="../reseq/cnv/segmented_coverage/1020540_Remission.segmented-coverage.tsv", circos="~/hdall/results/cnv/circos/1020540_dia.cnv.circos.tsv", stringsAsFactors=F)

# normalization factor (correct for different read count)
normal.chrs <- list(
		"243_dia" = c("chr2", "chr3", "chr16"),
		"243_rel" = c("chr1", "chr2", "chr3"),
		"933_dia" = c("chr1", "chr2", "chr3"), 
		"933_rel" = c("chr4", "chr9", "chr16"), 
		"944_dia" = c("chr1", "chr2", "chr3"), 
		"944_rel" = c("chr1", "chr2", "chr3"), 
		"1187_dia" = c("chr2", "chr3", "chr15"), 
		"1187_rel" = c("chr2", "chr3", "chr15"), 
		"1004564_dia" = c("chr1", "chr2", "chr3"),
		"1009302_dia" = c("chr1", "chr2", "chr3"), 
		"1009302_rel" = c("chr2", "chr3", "chr15"),
		"1010661_dia" = c("chr2", "chr3", "chr15"), 
		"1010781_dia" = c("chr2", "chr3", "chr5"),
		"1017005_rel" = c("chr2", "chr3", "chr16"),
		"1019357_dia" = c("chr1", "chr2", "chr3"),
		"1019357_rel" = c("chr1", "chr2", "chr3"), 
		"1019964_dia" = c("chr3", "chr11", "chr19"),
		"1020076_dia" = c("chr2", "chr3", "chr4"),
		"1020540_dia" = c("chr7", "chr8", "chr9", "chr10"),
		"1020540_rel" = c("chr1", "chr2", "chr3"),
		"1020583_dia" = c("chr2", "chr3", "chr16"),
		"1020583_rel" = c("chr1", "chr2", "chr3"),
		"1021087_dia" = c("chr1", "chr2", "chr3"),
		"1021392_dia" = c("chr1", "chr2", "chr3"),
		"1021392_rel" = c("chr1", "chr2", "chr3"),
		"1021631_dia" = c("chr2", "chr3", "chr16"),
		"1021631_rel" = c("chr2", "chr3", "chr16"),
		"1021865_rel" = c("chr2", "chr3", "chr16"),
		"1022914_dia" = c("chr13", "chr15", "chr16"),
		"1022914_rel" = c("chr13", "chr15", "chr16"),  # EXCLUDE
		"1023056_dia" = c("chr2", "chr3", "chr16"),
		"1023056_rel" = c("chr2", "chr3"), # EXCLUDE
		"1023338_dia" = c("chr2", "chr3", "chr16"),
		"1023392_dia" = c("chr3", "chr15", "chr16"),
		"1023392_rel" = c("chr1", "chr15", "chr16"),
		"1023545_rel" = c("chr2", "chr3", "chr15"),
		"1023616_dia" = c("chr2", "chr3", "chr16"),
		"1024518_dia" = c("chr2", "chr3", "chr16"),
		"1024518_rel" = c("chr2", "chr3", "chr15"),
		"1024543_dia" = c("chr2", "chr15"),
		"1024543_rel" = c("chr2", "chr15"),
		"1024589_dia" = c("chr2", "chr3", "chr16"),
		"1025108_dia" = c("chr2", "chr3", "chr16"),
		"1025108_rel" = c("chr3", "chr15", "chr16"),
		"1025409_dia" = c("chr2", "chr3", "chr16"),
		"1025409_rel" = c("chr2", "chr3", "chr16"),
		"1026233_dia" = c("chr2", "chr3", "chr16"),
		"1026662_dia" = c("chr2", "chr3", "chr15"),
		"AD15_rel" = c("chr1", "chr2", "chr3"),
		"B100_dia" = c("chr3", "chr7", "chr15"),
		"BL16_rel" = c("chr1", "chr2", "chr3"),  # EXCLUDE
		"BM18_rel" = c("chr2", "chr3"), # EXCLUDE
		"CA18_rel" = c("chr2", "chr3"), # EXCLUDE
		"DM1_rel" = c("chr2", "chr3", "chr16"),
		"EF7_dia" = c("chr3", "chr4", "chr15"),
		"FB11_dia" = c("chr2", "chr3", "chr15"),
		"FB11_rel" = c("chr2", "chr3", "chr15"),
		"FB14_dia" = c("chr2", "chr3", "chr15"),
		"FE1_rel" = c("chr2", "chr3", "chr15"),
		"FS1_rel" = c("chr2", "chr3"), # EXCLUDE
		"G44_dia" = c("chr3", "chr15", "chr16"),
		"GD1_rel" = c("chr2", "chr3", "chr15"),  # REINSPECT; cannot open pdf
		"GD18_rel" = c("chr2", "chr3", "chr15"),
		"G_dia" = c("chr2", "chr3", "chr15"),
		"G_rel" = c("chr2", "chr3"),
		"HD7_dia" = c("chr2", "chr3", "chr15"),
		"HJ15_rel" = c("chr2", "chr3", "chr15"),
		"HJA15_rel" = c("chr2", "chr3", "chr15"),
		"HL1_rel" = c("chr2", "chr3", "chr15"),
		"HS6_dia" = c("chr2", "chr3"),
		"HS6_rel" = c("chr2", "chr3", "chr15"),
		"KA17_rel" = c("chr2", "chr3", "chr15"),
		"KA14651_dia" = c("chr2", "chr3", "chr15"),
		"KA14651_rel" = c("chr2", "chr3", "chr15"),
		"KD20493_dia" = c("chr2", "chr3", "chr15"),
		"KD20493_rel" = c("chr2", "chr3", "chr15"),
		"K_dia" = c("chr2", "chr3", "chr15"),
		"K_rel" = c("chr2", "chr3", "chr16"), 
		"KE12025_dia" = c("chr2", "chr3"),
		"KE12025_rel" = c("chr2", "chr3", "chr15"),
		"KJ17_rel" = c("chr3", "chr15", "chr16"), # EXCLUDE
		"KL16_rel" = c("chr2", "chr3", "chr15"), 
		"LB17_rel" = c("chr2", "chr5", "chr7"), 
		"LM18_rel" = c("chr2", "chr3"), 
		"MB2_dia" = c("chr2", "chr3", "chr15"), 
		"MB2_rel" = c("chr2", "chr3", "chr15"), 
		"MJ1_rel" = c("chr2", "chr3", "chr15"), 
		"MJ16441_dia" = c("chr2", "chr3", "chr15"), 
		"MJ16441_rel" = c("chr2", "chr3", "chr15"), 
		"MV16_rel" = c("chr2", "chr3", "chr8"), 
		"NH17331_dia" = c("chr1", "chr2", "chr3"), 
		"NH17331_rel" = c("chr2", "chr3", "chr15"), 
		"NS18_rel" = c("chr2", "chr4"), # EXCLUDE
		"PC16_rel" = c("chr2", "chr3"), # EXCLUDE 
		"PJ13414_dia" = c("chr2", "chr3", "chr15"), 
		"PJ13414_rel" = c("chr2", "chr3", "chr15"), 
		"RD17412_dia" = c("chr1", "chr2", "chr3"), 
		"RD17412_rel" = c("chr1", "chr2", "chr3"), 
		"RS13466_dia" = c("chr2", "chr3", "chr15"), 
		"RS13466_rel" = c("chr2", "chr3", "chr15"), 
		"RT15_rel" = c("chr2", "chr3", "chr16"), 
		"RT16_rel" = c("chr2", "chr3", "chr16"), 
		"SJM16_rel" = c("chr2", "chr3"), # EXCLUDE 
		"SKR1_rel" = c("chr2", "chr3"), # EXCLUDE 
		"SL1_rel" = c("chr2", "chr3", "chr15"), 
		"SLM1_rel" = c("chr2", "chr3", "chr16"), 
		"ST14_rel" = c("chr2", "chr3", "chr15"), 
		"ST13892_dia" = c("chr2", "chr3", "chr15"), 
		"ST13892_rel" = c("chr2", "chr3", "chr15"), 
		"WA1_rel" = c("chr2", "chr3", "chr15"), 
		"ZA16211_dia" = c("chr1", "chr2", "chr3"), 
		"ZA16211_rel" = c("chr2", "chr7", "chr16"), 
		"ZE13_rel" = c("chr2", "chr3", "chr15"), 

		"39_dia" = c("chr2", "chr3"),
		"45_dia" = c("chr2", "chr8"),
		"49_dia" = c("chr2", "chr3", "chr16"),
		"54_dia" = c("chr2", "chr3"),
		"60_dia" = c("chr2", "chr3"),
		"73_dia" = c("chr2", "chr15"),
		"110_dia" = c("chr2", "chr3"),
		"111_dia" = c("chr2", "chr3", "chr15"),
		"134_dia" = c("chr2", "chr3", "chr12"),
		"143_dia" = c("chr2", "chr3", "chr15"),
		"147_dia" = c("chr2", "chr3", "chr15"),
		"199_dia" = c("chr2", "chr3", "chr15"),
		"331_dia" = c("chr2", "chr3", "chr16"),
		"350_dia" = c("chr2", "chr3", "chr15"),
		"380_dia" = c("chr2", "chr3", "chr15"),
		"409_dia" = c("chr2", "chr3", "chr15"),
		"442_dia" = c("chr2", "chr15", "chr16"),
		"461_dia" = c("chr2", "chr3", "chr16"),
		"466_dia" = c("chr9", "chr15"),
		"529_dia" = c("chr2", "chr3"),
		"530_dia" = c("chr2", "chr3", "chr16"),
		"591_dia" = c("chr2", "chr15"),
		"594_dia" = c("chr2", "chr15"),
		"602_dia" = c("chr2", "chr5"),
		"619_dia" = c("chr2", "chr3", "chr5"),
		"633_dia" = c("chr2", "chr15"),
		"634_dia" = c("chr2", "chr15"),
		"642_dia" = c("chr2", "chr5"),
		"646_dia" = c("chr2", "chr15"),
		"653_dia" = c("chr2", "chr3"),
		"666_dia" = c("chr2", "chr3"), # IKZF1 deletion
		"672_dia" = c("chr2", "chr5"),
		"681_dia" = c("chr2", "chr5"),
		"687_dia" = c("chr2", "chr7"),
		"697_dia" = c("chr2", "chr3"), # low blast
		"698_dia" = c("chr2", "chr3"), # hyperdiploid?
		"700_dia" = c("chr2", "chr3"),
		"709_dia" = c("chr2", "chr3"), # EXCLUDE
		"718_dia" = c("chr2", "chr3"),
		"724_dia" = c("chr2", "chr3"),
		"748_dia" = c("chr2", "chr3"),
		"754_dia" = c("chr2", "chr3"),
		"762_dia" = c("chr2", "chr15"),
		"776_dia" = c("chr2", "chr3"),
		"777_dia" = c("chr1", "chr2", "chr3"),
		"779_dia" = c("chr2", "chr3"),
		"782_dia" = c("chr2", "chr3"),
		"NRD_1_dia" = c("chr3", "chr8"),
		"NRD_2_dia" = c("chr2", "chr3"),
		"NRD_3_dia" = c("chr2", "chr3"),
		"NRD_4_dia" = c("chr2", "chr3"))

t <- read.delim(opt$tumor, header=F)
n <- read.delim(opt$normal, header=F)
m <- merge(n, t, by=c("V1", "V2"))
m$ratio <- log((m[,"V3.y"] / sum(m[,"V3.y"])) / (m[,"V3.x"] / sum(m[,"V3.x"])), 2)
m$ratio <- m$ratio - ifelse(!is.null(normal.chrs[[opt$patient]]), mean(m$ratio[which(m$V1 %in% normal.chrs[[opt$patient]] & is.finite(m$ratio))]), 0)

set.seed(25)
CNA.object <-CNA(genomdat = m[,"ratio"], chrom = m[,"V1"], maploc = m[,"V2"], data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
#segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2, undo.splits="sdundo", undo.SD=3)
segs <- segment(CNA.smoothed, alpha = 0.01, verbose=0, min.width=2)

pdf(opt$output, width=30, paper="A4r")
par(mfrow=c(5,5), mar=c(0.5,2,1.5,1))
for (c in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"))
{
	plot(m[m$V1==c,"V2"], m[m$V1==c,"ratio"], ylim=c(-1.5,1.5), main=c, cex=0.3, xaxt='n')
	abline(0, 0)
	for(i in which(segs$output$chrom==c)) {
		lines(c(segs$output$loc.start[i], segs$output$loc.end[i]),c(segs$output$seg.mean[i],segs$output$seg.mean[i]),type="l", col="orange", lty=1, lwd=3)
	}
}
dev.off()

# write circos output file
segs$output <- segs$output[complete.cases(segs$output),]
segs$output$color <- NA
segs$output$color[segs$output$seg.mean<0.15 & segs$output$seg.mean > -0.15] <- "fill_color=white"
segs$output$color[segs$output$seg.mean>=0.15 & segs$output$seg.mean<0.65] <- "fill_color=lblue"
segs$output$color[segs$output$seg.mean>=0.65] <- "fill_color=dblue"
segs$output$color[segs$output$seg.mean>-0.65 & segs$output$seg.mean<=-0.15] <- "fill_color=lred"
segs$output$color[segs$output$seg.mean<=-0.65] <- "fill_color=dred"
segs$output$seg.mean <- abs(segs$output$seg.mean)
segs$output$chrom <- gsub("chr", "hs", segs$output$chrom)
write.table(segs$output[!is.na(segs$output$seg.mean),c("chrom", "loc.start", "loc.end", "seg.mean", "color")], file=opt$circos, row.names=F, col.names=F, sep="\t", quote=F)
