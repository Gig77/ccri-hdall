library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)

rm(list=ls())

# assign colors to mutations
#cols <- c("KRAS:chr12:25378561:G>A" = brewer.pal(9, "YlOrRd")[2],
#		"KRAS:chr12:25378562:C>T" = brewer.pal(9, "YlOrRd")[3], 
#		"KRAS:chr12:25378562:C>A" = brewer.pal(9, "YlOrRd")[3], 
#		"KRAS:chr12:25378562:C>G" = brewer.pal(9, "YlOrRd")[3], 
#		"KRAS:chr12:25380275:T>A" = brewer.pal(9, "YlOrRd")[4],
#	    "KRAS:chr12:25380275:T>G" = brewer.pal(9, "YlOrRd")[4],
#		"KRAS:chr12:25398281:C>T" = brewer.pal(9, "YlOrRd")[5],
#		"KRAS:chr12:25398282:C>A" = brewer.pal(9, "YlOrRd")[5],
#		"KRAS:chr12:25398284:C>T" = brewer.pal(9, "YlOrRd")[6],
#		"KRAS:chr12:25398284:C>A" = brewer.pal(9, "YlOrRd")[6],
#		"KRAS:chr12:25398284:C>G" = brewer.pal(9, "YlOrRd")[6],
#		"KRAS:chr12:25398285:C>T" = brewer.pal(9, "YlOrRd")[7],
#		"KRAS:chr12:25398285:C>A" = brewer.pal(9, "YlOrRd")[7],
#		"NRAS:chr1:115256528:T>A" = brewer.pal(9, "Blues")[2],
#		"NRAS:chr1:115256528:T>G" = brewer.pal(9, "Blues")[2],
#		"NRAS:chr1:115256529:T>C" = brewer.pal(9, "Blues")[3],
#		"NRAS:chr1:115256530:G>T" = brewer.pal(9, "Blues")[4],
#		"NRAS:chr1:115258744:C>T" = brewer.pal(9, "Blues")[5],
#		"NRAS:chr1:115258744:C>A" = brewer.pal(9, "Blues")[5],
#		"NRAS:chr1:115258745:C>G" = brewer.pal(9, "Blues")[6],
#		"NRAS:chr1:115258745:C>A" = brewer.pal(9, "Blues")[6],
#		"NRAS:chr1:115258747:C>T" = brewer.pal(9, "Blues")[7],
#		"NRAS:chr1:115258747:C>G" = brewer.pal(9, "Blues")[7],
#		"NRAS:chr1:115258747:C>A" = brewer.pal(9, "Blues")[7],
#		"NRAS:chr1:115258748:C>A" = brewer.pal(9, "Blues")[8],
#		"NRAS:chr1:115258748:C>T" = brewer.pal(9, "Blues")[8],
#		"PTPN11:chr12:112888163:G>T" = brewer.pal(9, "RdPu")[2],
#		"PTPN11:chr12:112888163:G>C" = brewer.pal(9, "RdPu")[2],
#		"PTPN11:chr12:112888165:G>T" = brewer.pal(9, "RdPu")[3],
#		"PTPN11:chr12:112888166:A>G" = brewer.pal(9, "RdPu")[3],
#		"PTPN11:chr12:112888189:G>A" = brewer.pal(9, "RdPu")[4],
#		"PTPN11:chr12:112888198:G>A" = brewer.pal(9, "RdPu")[4],
#		"PTPN11:chr12:112888198:G>T" = brewer.pal(9, "RdPu")[4],
#		"PTPN11:chr12:112888199:C>A" = brewer.pal(9, "RdPu")[4],
#		"PTPN11:chr12:112888199:C>T" = brewer.pal(9, "RdPu")[4],
#		"PTPN11:chr12:112888210:G>A" = brewer.pal(9, "RdPu")[5],
#		"PTPN11:chr12:112888211:A>G" = brewer.pal(9, "RdPu")[5],
#		"PTPN11:chr12:112888211:A>T" = brewer.pal(9, "RdPu")[5],
#		"PTPN11:chr12:112915524:A>G" = brewer.pal(9, "RdPu")[6],
#		"PTPN11:chr12:112926884:T>C" = brewer.pal(9, "RdPu")[7],
#		"PTPN11:chr12:112926884:T>A" = brewer.pal(9, "RdPu")[7],
#		"PTPN11:chr12:112926884:T>G" = brewer.pal(9, "RdPu")[7],
#		"PTPN11:chr12:112926888:G>T" = brewer.pal(9, "RdPu")[8],
#		"PTPN11:chr12:112926888:G>C" = brewer.pal(9, "RdPu")[8],
#		"PTPN11:chr12:112926909:A>T" = brewer.pal(9, "RdPu")[9],
#		"FLT3:chr13:28592634:CATG>C" = "black",
#		"FLT3:chr13:28592640:A>T" = "black",
#		"FLT3:chr13:28592640:A>T" = "black",
#		"FLT3:chr13:28592641:T>A" = "black",
#		"FLT3:chr13:28592642:C>A" = "black",
#		"FLT3:chr13:28592642:C>G" = "black",
#		"FLT3:chr13:28597522:C>T" = "black",
#		"FLT3:chr13:28592621:A>G" = "black",
#		"FLT3:chr13:28592623:T>G" = "black",
#		"FLT3:chr13:28608341:T>C" = "black",
#		"FLT3:chr13:28608312:T>G" = "black",
#		"FLT3:chr13:28592623:T>G" = "black",   
#		"FLT3:chr13:28592617:A>T" = "black",
#		"FLT3:chr13:28602329:G>A" = "black",
#		"FLT3:chr13:28602340:G>T" = "black",  
#		"FLT3:chr13:28602378:T>C" = "black",  
#		"FLT3:chr13:28608098:T>A" = "black",   
#		"FLT3:chr13:28592641:T>A" = "black",
#		"FLT3:chr13:28592623:T>G" = "black",   
#		"FLT3:chr13:28602329:G>A" = "black",
#		"FLT3:chr13:28608285:A>C" = "black",
#		"FLT3:chr13:28608341:T>C" = "black",
#		"FLT3:chr13:28589780:A>C" = "black",
#		"FLT3:chr13:28608329:A>T" = "black",
#		"FLT3:chr13:28610138:G>A" = "black",
#		"FLT3:chr13:28602380:T>C" = "black",
#		"FLT3:chr13:28608275:A>G" = "black",
#		"FLT3:chr13:28608320:A>C" = "black",
#		"FLT3:chr13:28608449:T>C" = "black",
#		"FLT3:chr13:28622527:A>G" = "black",
#		"FLT3:chr13:28626716:C>T" = "black",
#		"FLT3:chr13:28592717:T>C" = "black",
#		"FLT3:chr13:28608291:A>C" = "black")

patients.matched = c("314", "1021247", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "A", "B", "C", "D", "X", "Y", "1009302", "1019357", "1020540", "1020583", "1021392", "1021631", "1022914", "1023056", "1023392", "1024518", "1024543", "1025108", "1025409", "1187", "FB11", "G", "HS6", "K", "MB2", "243", "933", "944", "KA14651", "KD20493", "KE12025", "MJ16441", "NH17331", "PJ13414", "RD17412", "RS13466", "ST13892", "ZA16211")

cols <- c(		"KRAS" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25378561:G>A" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25378562:C>T" = brewer.pal(8, "Accent")[1], 
				"KRAS:chr12:25378562:C>A" = brewer.pal(8, "Accent")[1], 
				"KRAS:chr12:25378562:C>G" = brewer.pal(8, "Accent")[1], 
				"KRAS:chr12:25380275:T>A" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25380275:T>G" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398281:C>T" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398282:C>A" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398284:C>T" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398284:C>A" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398284:C>G" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398285:C>T" = brewer.pal(8, "Accent")[1],
				"KRAS:chr12:25398285:C>A" = brewer.pal(8, "Accent")[1],
				
				"NRAS" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115256528:T>A" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115256528:T>G" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115256529:T>C" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115256530:G>T" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258744:C>T" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258744:C>A" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258745:C>G" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258745:C>A" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258747:C>T" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258747:C>G" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258747:C>A" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258748:C>A" = brewer.pal(8, "Accent")[2],
				"NRAS:chr1:115258748:C>T" = brewer.pal(8, "Accent")[2],
				
				"PTPN11" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888163:G>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888163:G>C" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888165:G>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888166:A>G" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888189:G>A" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888198:G>A" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888198:G>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888199:C>A" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888199:C>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888210:G>A" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888211:A>G" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112888211:A>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112915524:A>G" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926884:T>C" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926884:T>A" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926884:T>G" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926888:G>T" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926888:G>C" = brewer.pal(8, "Accent")[5],
				"PTPN11:chr12:112926909:A>T" = brewer.pal(8, "Accent")[5],
				
				"FLT3" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592634:CATG>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592640:A>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592640:A>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592641:T>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592642:C>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592642:C>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28597522:C>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592621:A>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592623:T>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608341:T>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608312:T>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592623:T>G" = brewer.pal(8, "Accent")[7],   
				"FLT3:chr13:28592617:A>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28602329:G>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28602340:G>T" = brewer.pal(8, "Accent")[7],  
				"FLT3:chr13:28602378:T>C" = brewer.pal(8, "Accent")[7],  
				"FLT3:chr13:28608098:T>A" = brewer.pal(8, "Accent")[7],   
				"FLT3:chr13:28592641:T>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592623:T>G" = brewer.pal(8, "Accent")[7],   
				"FLT3:chr13:28602329:G>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608285:A>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608341:T>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28589780:A>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608329:A>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28610138:G>A" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28602380:T>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608275:A>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608320:A>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608449:T>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28622527:A>G" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28626716:C>T" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28592717:T>C" = brewer.pal(8, "Accent")[7],
				"FLT3:chr13:28608291:A>C" = brewer.pal(8, "Accent")[7],
				
				"diagnosis" = brewer.pal(8, "Accent")[1],
				"relapse" = brewer.pal(8, "Accent")[2])
				
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
c <- data.frame(patient=factor(c$patient_id[c$patient_id %in% m$patient], levels=levels(m$patient)), 
		        blasts.dia = suppressWarnings(as.numeric(c$blasts_dia[c$patient_id %in% m$patient])), 
				blasts.rel = suppressWarnings(as.numeric(c$blasts_rel[c$patient_id %in% m$patient])))
m <- merge(m, c, all.x=T)

m$frequency.norm <- NA
m$frequency.norm[m$sample=="diagnosis"] <- m$frequency[m$sample=="diagnosis"]/m$blasts.dia[m$sample=="diagnosis"]*100
m$frequency.norm[m$sample=="relapse"] <- m$frequency[m$sample=="relapse"]/m$blasts.rel[m$sample=="relapse"]*100
m$frequency.norm[is.na(m$frequency.norm)] <- m$frequency[is.na(m$frequency.norm)]
m$patient.label <- as.factor(ifelse(m$patient %in% patients.matched, paste0(m$patient,"*"), as.character(m$patient)))

# subset samples
m.relapsing <- m[m$cohort=="relapsing",]
write.table(m[m$cohort=="relapsing",], "~/hdall/results/ras-heterogeneity/ras-pathway-heterogeneity.mutations.relapsing-cohort.tsv", row.names=F, col.names=T, sep="\t", quote=F)

m.relapsing.dia <- m.relapsing[m.relapsing$sample=="diagnosis",]
m.relapsing.dia <- ddply(m.relapsing.dia[order(m.relapsing.dia$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.relapsing.rel <- m.relapsing[m.relapsing$sample=="relapse",]
m.relapsing.rel <- ddply(m.relapsing.rel[order(m.relapsing.rel$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.relapsing.cons <- merge(m.relapsing.dia, m.relapsing.rel, by=c("patient", "mut"))[,c("patient", "mut", "frequency.norm.x", "frequency.norm.y")]
m.relapsing.cons <- paste0(m.relapsing.cons$patient, ":", m.relapsing.cons$mut)
m.nonrelapsing <- m[m$cohort=="non-relapsing" & m$sample=="diagnosis",]

stop("OK")

pdf("~/hdall/results/ras-heterogeneity/ras-barplot.pdf", width=18, height=13)

#---
# RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient.label, data=m.relapsing.dia, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient.label", "frequency.norm")]
m.relapsing.dia$patient.label <- factor(m.relapsing.dia$patient.label, sorted$patient.label)

# find heterogenious samples
het.dia <- aggregate(frequency.norm~patient, data=m.relapsing.dia, FUN=length)
het.dia.patient <- het.dia[het.dia$frequency.norm>1,"patient"]
#m.relapsing.dia$group <- ifelse(m.relapsing.dia$patient %in% het.dia.patient, sprintf("heterogenious (%.0f%%)", length(het.dia.patient)/length(levels(m.relapsing.dia$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.dia.patient)/length(levels(m.relapsing.dia$patient))*100))
m.relapsing.dia$group <- ifelse(m.relapsing.dia$patient %in% het.dia.patient, "heterogeneous", "homogeneous")
m.relapsing.dia$label <- ifelse(paste0(m.relapsing.dia$patient, ":", m.relapsing.dia$mut) %in% m.relapsing.cons, "symbol(\"\\267\")", NA)

plot.dia <- ggplot(data=m.relapsing.dia, aes(x=patient.label, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(.~group, scale="free_x", space = "free_x") +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
		scale_fill_manual(values = cols, breaks=c("KRAS:chr12:25378561:G>A", "NRAS:chr1:115256528:T>A", "PTPN11:chr12:112888163:G>T", "FLT3:chr13:28589780:A>C"), labels=c("KRAS", "NRAS", "PTPN11", "FLT3")) + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.title = element_blank(), legend.text=element_text(size=13), legend.key.size=unit(0.7, "cm"), legend.justification=c(1,1), legend.position=c(1,1), legend.background=element_rect(colour="black"),
				strip.text.x = element_text(size=15), 
				panel.grid.major.x = element_blank(),
				panel.grid.minor.y = element_blank(),
				plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(0.5,1,0,1), "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Allelic frequency diagnosis") 
#		+ ggtitle("Ras pathway mutations at diagnosis")

#---
# RELAPSING, REL
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient", "frequency.norm")]
m.relapsing.rel$patient <- factor(m.relapsing.rel$patient, sorted$patient)

# find heterogenious samples
het.rel <- aggregate(frequency.norm~patient, data=m.relapsing.rel, FUN=length)
het.rel.patient <- het.rel[het.rel$frequency.norm>1,"patient"]
#m.relapsing.rel$group <- ifelse(m.relapsing.rel$patient %in% het.rel.patient, sprintf("heterogenious (%.0f%%)", length(het.rel.patient)/length(levels(m.relapsing.rel$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.rel.patient)/length(levels(m.relapsing.rel$patient))*100))
m.relapsing.rel$group <- ifelse(m.relapsing.rel$patient %in% het.rel.patient, "heterogeneous", "homogeneous")
m.relapsing.rel$label <- ifelse(paste0(m.relapsing.rel$patient, ":", m.relapsing.rel$mut) %in% m.relapsing.cons, "symbol(\"\\267\")", NA)

plot.rel <- ggplot(data=m.relapsing.rel, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
				scale_fill_manual(values = cols, name="Mutation", guide=FALSE) + 
				scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.x = element_text(size=20, vjust=0.2), axis.title.y = element_text(size=18, vjust=0.1), 
						legend.key.size=unit(0.35, "cm"),
						strip.text.x = element_text(size=15), 
						panel.grid.major.x = element_blank(),
						panel.grid.minor.y = element_blank(),
						plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(1,1,0.3,1), "cm")) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Allelic frequency relapse") 
#				+ ggtitle("Ras pathway mutations at relapse")
		
grid.arrange(plot.dia, plot.rel, nrow=2)

#---
# P-value
#---
sig <- fisher.test(matrix(c(sum(het.dia$frequency.norm>1), sum(het.dia$frequency.norm==1), sum(het.rel$frequency.norm>1), sum(het.rel$frequency.norm==1)), nrow=2, byrow=T))
print(sprintf("P-value (Fisher's exact test): %g", sig$p.value))

dev.off()

stop("OK")

#---
# RELAPSING, MATCHED, DIA+REL
#---

# reorder patients by cumulative frequency
m.relapsing.combined <- m.relapsing.dia[m.relapsing.dia$patient %in% patients.matched,]
m.relapsing.combined <- rbind(m.relapsing.combined, m.relapsing.rel[m.relapsing.rel$patient %in% patients.matched,])

# reorder patients by group and cumulative frequency
sorted <- aggregate(frequency.norm~patient+sample+group, data=m.relapsing.combined, FUN=sum)
sorted$sample <- as.character(sorted$sample)
sorted$sample[sorted$sample=="relapse"] <- "0relapse"
sorted$group[sorted$group=="homogeneous"] <- "0homogeneous"
sorted <- sorted[order(sorted$sample, sorted$group, sorted$frequency.norm, decreasing=T),]
m.relapsing.combined$patient <- factor(m.relapsing.combined$patient, unique(as.character(sorted$patient)))

pdf("~/hdall/results/ras-heterogeneity/ras-barplot-matched.pdf", width=14, height=9)
plot.dia <- ggplot(data=m.relapsing.combined, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(sample~.) +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
		scale_fill_manual(values = cols, breaks=c("KRAS:chr12:25378561:G>A", "NRAS:chr1:115256528:T>A", "PTPN11:chr12:112888163:G>T", "FLT3:chr13:28589780:A>C"), labels=c("KRAS", "NRAS", "PTPN11", "FLT3")) + 
		scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.title = element_blank(), legend.text=element_text(size=13), legend.key.size=unit(0.7, "cm"), legend.justification=c(1,1), legend.position=c(1,1), legend.background=element_rect(colour="black"),
				strip.text.x = element_text(size=15), 
				panel.grid.major.x = element_blank(),
				panel.grid.minor.y = element_blank(),
				plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(0.5,1,0,1), "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Allelic frequency") 
#		+ ggtitle("Ras pathway mutations at diagnosis")
print(plot.dia)
dev.off()

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
m.nonrelapsing$group <- ifelse(m.nonrelapsing$patient %in% het, "heterogeneous", "homogeneous")

pdf("~/hdall/results/ras-heterogeneity/ras-barplot-nonrelapsing.pdf", width=14, height=9)
print(ggplot(data=m.nonrelapsing, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) +
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				scale_fill_manual(values = cols, breaks=c("KRAS:chr12:25378561:G>A", "NRAS:chr1:115256528:T>G", "PTPN11:chr12:112888163:G>T", "FLT3:chr13:28592642:C>A"), labels=c("KRAS", "NRAS", "PTPN11", "FLT3"), name="") + 
#				scale_fill_manual(values = cols, name="Mutation") + 
				scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks=c(0.25, 0.5, 0.75, 1)) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1), 
						legend.key.size=unit(0.35, "cm"),
						panel.grid.major.x = element_blank(),
						panel.grid.minor.y = element_blank()) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Allelic frequency") +
				ggtitle("RAS pathway mutations at diagnosis (non-relapsing cohort)"))
dev.off()

