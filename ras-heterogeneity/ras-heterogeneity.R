rm(list=ls())

hotspots <- data.frame(gene=c("KRAS", "KRAS", "KRAS", "KRAS", "KRAS", "KRAS", "KRAS", "NRAS", "NRAS", "NRAS", "NRAS", "NRAS", "NRAS", "NRAS", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "PTPN11", "FLT3", "FLT3", "FLT3", "FLT3"), 
		  			   chr=c("chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr12", "chr13", "chr13", "chr13", "chr13"),
					   pos=c(25398281, 25398282, 25398284, 25398285, 25378561, 25378562, 25380275, 115256528, 115256529, 115256530, 115258744, 115258745, 115258747, 115258748, 112888162, 112888163, 112888165, 112888166, 112888189, 112888197, 112888198, 112888199, 112888210, 112888211, 112915523, 112915524, 112926884, 112926885, 112926888, 112926908, 112926909, 28592640, 28592641, 28592642, 28592621),
					   stringsAsFactors=F)
#hotspots <- data.frame(gene=c("PTPN11"), chr=c("chr12"), pos=c(112888189), stringsAsFactors=F)			   
mutations <- data.frame(patient=character(0), cohort=character(0), sample=character(0), gene=character(0), chr=character(0), pos=numeric(0), ref=character(0), alt=character(0), alt.reads=numeric(0), tot.reads=numeric(0), frequency=numeric(0), stringsAsFactors=F)

patients.matched <- c("KA14651", "1009302", "1019357", "1020540", "1020583", "1021247", "1021392", "1021631", "1022914", "1023056", "1023392", "1024518", "1024543", "1025108", "1025409", "1187", "314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "A", "B", "C", "D", "FB11", "G", "HS6", "K", "MB2", "X", "243", "933", "944", "KD20493", "KE12025", "MJ16441", "NH17331", "PJ13414", "RD17412", "RS13466", "ST13892", "ZA16211", "Y")
patients.relonly <- c("1017005", "1021865", "1023545", "AD15", "BL16", "BM18", "CA18", "DM1", "FE1", "FS1", "GD18", "GD1", "HJ15", "HJA15", "HL1", "KA17", "KJ17", "KL16", "LB17", "LM18", "MJ1", "ML10", "MV16", "NS18", "PC16", "RT15", "RT16", "SJM16", "SKR1", "SL1", "SLM1", "ST14", "WA1", "ZE13")
patients.diaonly <- c("1004564", "1010661", "1010781", "1019964", "1020076", "1021087", "1023338", "1023616", "1024589", "1026233", "1026662", "B100", "EF7", "FB14", "G44", "HD7")
patients.dianonrel <- c("331", "380", "442", "350", "461", "466", "529", "591", "602", "619", "633", "634", "642", "653", "666", "672", "697", "698", "700", "709", "724", "762", "776", "777", "779", "782", "409", "NRD_1", "73", "NRD_2", "NRD_3", "60", "594", "687", "748", "754", "646", "530", "718", "681", "39", "49", "45", "54", "110", "111", "134", "143", "147", "199", "NRD_4")

check_base <- function(p, cohort, sample, gene, chr, pos, ref, base, qual, mutations) {
	for (b in c("A", "T", "G", "C")) {
		#print(paste0("ALT", ": ", sum(base==tolower(b) | base==toupper(b)), " ", mean(qual[base==tolower(b) | base==toupper(b)])))
		#print(paste0("REF", ": ", sum(base=="," | base=="."), " ", mean(qual[base=="," | base=="."])))
		if (sum(base==tolower(b)) >= 2 && sum(base==toupper(b)) >= 2 && mean(qual[base==tolower(b)]) >= 22 && mean(qual[base==toupper(b)]) >= 22)
		{
			num.alt <- sum(toupper(base)==b)
			tot <- length(base)
			freq <- num.alt/tot
			print(paste0(p, "\t", cohort, "\t", sample, "\t", gene, "\t", chr, ":", pos, "\t", ref, "->", b, "\tALT:", num.alt, "\tTOT:", tot, "\tFREQ:", freq))
			mutations[nrow(mutations)+1,] <- c(p, cohort, sample, gene, chr, pos, ref, b, num.alt, tot, freq)
		}	
	}
	
	return(mutations)
}

process_bam <- function(p, cohort, sample, bam, gene, chr, pos, mutations) {
	cmd <- paste0("samtools view -bh ", bam, " ", chr, ":", pos, "-", pos+1, " | samtools mpileup -q 1 -f /data/christian/generic/data/current/hg19/ucsc.hg19.fasta - 2>/dev/null | grep ", pos)
	#print(cmd)
	r <- system(cmd, intern=T)
	#print(r)
	
	ref <- strsplit(r, "\t")[[1]][3]
	base <- strsplit(r, "\t")[[1]][5]
	base <- gsub("\\^.", "", base, perl=T)
	base <- gsub("\\$", "", base, perl=T)
	base <- strsplit(base, "")[[1]]
	#print(base)
	qual <- as.numeric(charToRaw(strsplit(r, "\t")[[1]][6]))-33
	#print(qual)
	
	mutations <- check_base(p, cohort, sample, gene, chr, pos, ref, base, qual, mutations)
	
	return(mutations)
}

for (i in 1:nrow(hotspots)) {
	for (p in patients.dianonrel) {
		mutations <- process_bam(p, "non-relapsing", "diagnosis", paste0("~/hdall/data/reseq/bam/", p, "_Diagnosis.duplicate_marked.realigned.recalibrated.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
	for (p in c(patients.matched, patients.diaonly)) {
		mutations <- process_bam(p, "relapsing", "diagnosis", paste0("~/hdall/data/reseq/bam/", p, "_Diagnosis.duplicate_marked.realigned.recalibrated.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
	for (p in c(patients.matched, patients.relonly)) {
		mutations <- process_bam(p, "relapsing", "relapse", paste0("~/hdall/data/reseq/bam/", p, "_Relapse.duplicate_marked.realigned.recalibrated.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
	}
#	for (p in c("ZE13")) {
#		mutations <- process_bam(p, "relapsing", "relapse", paste0("~/hdall/data/reseq/bam/", p, "_Relapse.duplicate_marked.realigned.recalibrated.bam"), hotspots[i, "gene"], hotspots[i, "chr"], hotspots[i, "pos"], mutations)
#	}
}

write.table(mutations, file="~/hdall/results/ras-heterogeneity/ras.hotspots.tsv", col.names=T, row.names=F, sep="\t", quote=F)


