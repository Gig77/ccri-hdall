library(RColorBrewer)
library(ggplot2)
library(grid)

rm(list=ls())

# assign colors to mutations
cols <- c("KRAS:chr12:25378561:G>A" = brewer.pal(9, "YlOrRd")[2],
		"KRAS:chr12:25378562:C>T" = brewer.pal(9, "YlOrRd")[3], 
		"KRAS:chr12:25378562:C>A" = brewer.pal(9, "YlOrRd")[3], 
		"KRAS:chr12:25378562:C>G" = brewer.pal(9, "YlOrRd")[3], 
		"KRAS:chr12:25380275:T>A" = brewer.pal(9, "YlOrRd")[4],
	    "KRAS:chr12:25380275:T>G" = brewer.pal(9, "YlOrRd")[4],
		"KRAS:chr12:25398281:C>T" = brewer.pal(9, "YlOrRd")[5],
		"KRAS:chr12:25398282:C>A" = brewer.pal(9, "YlOrRd")[5],
		"KRAS:chr12:25398284:C>T" = brewer.pal(9, "YlOrRd")[6],
		"KRAS:chr12:25398284:C>A" = brewer.pal(9, "YlOrRd")[6],
		"KRAS:chr12:25398284:C>G" = brewer.pal(9, "YlOrRd")[6],
		"KRAS:chr12:25398285:C>T" = brewer.pal(9, "YlOrRd")[7],
		"KRAS:chr12:25398285:C>A" = brewer.pal(9, "YlOrRd")[7],
		"NRAS:chr1:115256528:T>A" = brewer.pal(9, "Blues")[2],
		"NRAS:chr1:115256528:T>G" = brewer.pal(9, "Blues")[2],
		"NRAS:chr1:115256529:T>C" = brewer.pal(9, "Blues")[3],
		"NRAS:chr1:115256530:G>T" = brewer.pal(9, "Blues")[4],
		"NRAS:chr1:115258744:C>T" = brewer.pal(9, "Blues")[5],
		"NRAS:chr1:115258744:C>A" = brewer.pal(9, "Blues")[5],
		"NRAS:chr1:115258745:C>G" = brewer.pal(9, "Blues")[6],
		"NRAS:chr1:115258745:C>A" = brewer.pal(9, "Blues")[6],
		"NRAS:chr1:115258747:C>T" = brewer.pal(9, "Blues")[7],
		"NRAS:chr1:115258747:C>G" = brewer.pal(9, "Blues")[7],
		"NRAS:chr1:115258747:C>A" = brewer.pal(9, "Blues")[7],
		"NRAS:chr1:115258748:C>A" = brewer.pal(9, "Blues")[8],
		"NRAS:chr1:115258748:C>T" = brewer.pal(9, "Blues")[8],
		"PTPN11:chr12:112888163:G>T" = brewer.pal(9, "RdPu")[2],
		"PTPN11:chr12:112888163:G>C" = brewer.pal(9, "RdPu")[2],
		"PTPN11:chr12:112888165:G>T" = brewer.pal(9, "RdPu")[3],
		"PTPN11:chr12:112888166:A>G" = brewer.pal(9, "RdPu")[3],
		"PTPN11:chr12:112888189:G>A" = brewer.pal(9, "RdPu")[4],
		"PTPN11:chr12:112888198:G>A" = brewer.pal(9, "RdPu")[4],
		"PTPN11:chr12:112888198:G>T" = brewer.pal(9, "RdPu")[4],
		"PTPN11:chr12:112888199:C>A" = brewer.pal(9, "RdPu")[4],
		"PTPN11:chr12:112888199:C>T" = brewer.pal(9, "RdPu")[4],
		"PTPN11:chr12:112888210:G>A" = brewer.pal(9, "RdPu")[5],
		"PTPN11:chr12:112888211:A>G" = brewer.pal(9, "RdPu")[5],
		"PTPN11:chr12:112888211:A>T" = brewer.pal(9, "RdPu")[5],
		"PTPN11:chr12:112915524:A>G" = brewer.pal(9, "RdPu")[6],
		"PTPN11:chr12:112926884:T>C" = brewer.pal(9, "RdPu")[7],
		"PTPN11:chr12:112926884:T>A" = brewer.pal(9, "RdPu")[7],
		"PTPN11:chr12:112926884:T>G" = brewer.pal(9, "RdPu")[7],
		"PTPN11:chr12:112926888:G>T" = brewer.pal(9, "RdPu")[8],
		"PTPN11:chr12:112926888:G>C" = brewer.pal(9, "RdPu")[8],
		"PTPN11:chr12:112926909:A>T" = brewer.pal(9, "RdPu")[9],
		"FLT3:chr13:28592634:CATG>C" = "black",
		"FLT3:chr13:28592640:A>T" = "black",
		"FLT3:chr13:28592640:A>T" = "black",
		"FLT3:chr13:28592641:T>A" = "black",
		"FLT3:chr13:28592642:C>A" = "black",
		"FLT3:chr13:28592642:C>G" = "black",
		"FLT3:chr13:28597522:C>T" = "black",
		"FLT3:chr13:28592621:A>G" = "black",
		"FLT3:chr13:28592623:T>G" = "black",
		"FLT3:chr13:28608341:T>C" = "black",
		"FLT3:chr13:28608312:T>G" = "black",
		"FLT3:chr13:28592623:T>G" = "black",   
		"FLT3:chr13:28592617:A>T" = "black",
		"FLT3:chr13:28602329:G>A" = "black",
		"FLT3:chr13:28602340:G>T" = "black",  
		"FLT3:chr13:28602378:T>C" = "black",  
		"FLT3:chr13:28608098:T>A" = "black",   
		"FLT3:chr13:28592641:T>A" = "black",
		"FLT3:chr13:28592623:T>G" = "black",   
		"FLT3:chr13:28602329:G>A" = "black",
		"FLT3:chr13:28608285:A>C" = "black",
		"FLT3:chr13:28608341:T>C" = "black",
		"FLT3:chr13:28589780:A>C" = "black",
		"FLT3:chr13:28608329:A>T" = "black",
		"FLT3:chr13:28610138:G>A" = "black",
		"FLT3:chr13:28602380:T>C" = "black",
		"FLT3:chr13:28608275:A>G" = "black",
		"FLT3:chr13:28608320:A>C" = "black",
		"FLT3:chr13:28608449:T>C" = "black",
		"FLT3:chr13:28622527:A>G" = "black",
		"FLT3:chr13:28626716:C>T" = "black",
		"FLT3:chr13:28592717:T>C" = "black",
		"FLT3:chr13:28608291:A>C" = "black")

# get KRAS, NRAS, PTPN11 mutations from hotspot mutation caller
m <- read.delim("~/hdall/results/ras-heterogeneity/ras.hotspots.tsv")
m <- m[m$gene != "FLT3",]

# get FLT3 mutations from MuTect
reseq.relapsing <- read.delim("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", stringsAsFactors=F)

flt3.relapsing.dia <- reseq.relapsing[reseq.relapsing$status!="REJECT" & reseq.relapsing$non_silent==T & reseq.relapsing$sample=="rem_dia" & reseq.relapsing$gene=="FLT3", c("patient", "sample", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
flt3.relapsing.dia$sample <- "diagnosis"
flt3.relapsing.dia <- data.frame(patient=flt3.relapsing.dia[,1], cohort="relapsing", flt3.relapsing.dia[,2:ncol(flt3.relapsing.dia)])
flt3.relapsing.dia$patient <- as.character(flt3.relapsing.dia$patient)
flt3.relapsing.dia$cohort <- as.character(flt3.relapsing.dia$cohort)
names(flt3.relapsing.dia) <- names(m)
m <- rbind(m, flt3.relapsing.dia)

flt3.relapsing.rel <- reseq.relapsing[reseq.relapsing$status!="REJECT" & reseq.relapsing$non_silent==T & reseq.relapsing$sample=="rem_rel" & reseq.relapsing$gene=="FLT3", c("patient", "sample", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
flt3.relapsing.rel$sample <- "relapse"
flt3.relapsing.rel <- data.frame(patient=flt3.relapsing.rel[,1], cohort="relapsing", flt3.relapsing.rel[,2:ncol(flt3.relapsing.rel)])
flt3.relapsing.rel$patient <- as.character(flt3.relapsing.rel$patient)
flt3.relapsing.rel$cohort <- as.character(flt3.relapsing.rel$cohort)
names(flt3.relapsing.rel) <- names(m)
m <- rbind(m, flt3.relapsing.rel)

reseq.nonrelapsing <- read.delim("~/hdall/results/reseq/filtered-variants.reseq.nonrel.tsv", stringsAsFactors=F)
reseq.nonrelapsing$freq_leu <- reseq.nonrelapsing$freq_leu / 100 
flt3.nonrelapsing <- reseq.nonrelapsing[reseq.nonrelapsing$status!="REJECT" & reseq.nonrelapsing$non_silent==T & reseq.nonrelapsing$gene=="FLT3", c("patient", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
flt3.nonrelapsing <- data.frame(patient=flt3.nonrelapsing[,1], cohort="non-relapsing", sample="diagnosis", flt3.nonrelapsing[,2:ncol(flt3.nonrelapsing)])
flt3.nonrelapsing$patient <- as.character(flt3.nonrelapsing$patient)
flt3.nonrelapsing$cohort <- as.character(flt3.nonrelapsing$cohort)
flt3.nonrelapsing$sample <- as.character(flt3.nonrelapsing$sample)
names(flt3.nonrelapsing) <- names(m)
m <- rbind(m, flt3.nonrelapsing)

# add columns
m$mut <- paste0(m$gene, ":", m$chr, ":", m$pos, ":", m$ref, ">", m$alt)
m$mut.short <- paste0(m$ref, ">", m$alt)

# normalize by blast counts
c <- read.delim("~/hdall/results/clinical/clinical_data.tsv", , na.strings=c("", "NA", "n/a", "n/d", " ", "early (CNS)"), stringsAsFactors=F)
c <- data.frame(patient=factor(c$patient_id, levels=levels(m$patient)), blasts.dia = as.numeric(c$blasts_dia), blasts.rel = as.numeric(c$blasts_rel))
m <- merge(m, c, all.x=T)

m$frequency.norm <- NA
m$frequency.norm[m$sample=="diagnosis"] <- m$frequency[m$sample=="diagnosis"]/m$blasts.dia[m$sample=="diagnosis"]*100
m$frequency.norm[m$sample=="relapse"] <- m$frequency[m$sample=="relapse"]/m$blasts.rel[m$sample=="relapse"]*100

# subset samples
m.relapsing.dia <- m[m$cohort=="relapsing" & m$sample=="diagnosis",]
m.relapsing.rel <- m[m$cohort=="relapsing" & m$sample=="relapse",]
m.nonrelapsing <- m[m$cohort=="non-relapsing" & m$sample=="diagnosis",]

pdf("~/hdall/results/ras-heterogeneity/ras-barplot.pdf", width=15)

#---
# RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.relapsing.dia, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.relapsing.dia$patient <- factor(m.relapsing.dia$patient, sorted$patient)

# find heterogenious samples
het <- aggregate(frequency.norm~patient, data=m.relapsing.dia, FUN=length)
het <- het[het$frequency.norm>1,"patient"]
m.relapsing.dia$group <- ifelse(m.relapsing.dia$patient %in% het, "heterogenious", "homogenious")

print(ggplot(data=m.relapsing.dia, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(.~group, scale="free_x", space = "free_x") +
		geom_bar(stat="identity", width=0.9, colour="white") +
		scale_fill_manual(values = cols, name="Mutation") + 
		scale_y_continuous(expand = c(0,0), limits = c(0,max(sorted$frequency.norm+0.1)), breaks=seq(0, 1, 0.1)) + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size=unit(0.35, "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Allelic frequency (normalized)") +
		ggtitle("RAS pathway mutations at diagnosis (relapsing cohort)"))

#---
# RELAPSING, REL
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.relapsing.rel$patient <- factor(m.relapsing.rel$patient, sorted$patient)

# find heterogenious samples
het <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=length)
het <- het[het$frequency.norm>1,"patient"]
m.relapsing.rel$group <- ifelse(m.relapsing.rel$patient %in% het, "heterogenious", "homogenious")

print(ggplot(data=m.relapsing.rel, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				scale_fill_manual(values = cols, name="Mutation") + 
				scale_y_continuous(expand = c(0,0), limits = c(0,max(sorted$frequency.norm+0.1)), breaks=seq(0, 1, 0.1)) + 
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size=unit(0.35, "cm")) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Allelic frequency (normalized)") +
				ggtitle("RAS pathway mutations at relapse (relapsing cohort)"))

#---
# NON-RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.nonrelapsing, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.nonrelapsing$patient <- factor(m.nonrelapsing$patient, sorted$patient)

# find heterogenious samples
het <- aggregate(frequency.norm~patient, data=m.nonrelapsing, FUN=length)
het <- het[het$frequency.norm>1,"patient"]
m.nonrelapsing$group <- ifelse(m.nonrelapsing$patient %in% het, "heterogenious", "homogenious")

print(ggplot(data=m.nonrelapsing, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) +
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				scale_fill_manual(values = cols, name="Mutation") + 
				scale_y_continuous(expand = c(0,0), limits = c(0,max(sorted$frequency.norm+0.1)), breaks=seq(0, 1, 0.1)) + 
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size=unit(0.35, "cm")) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Allelic frequency (normalized)") +
				ggtitle("RAS pathway mutations at diagnosis (non-relapsing cohort)"))


dev.off()

