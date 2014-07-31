library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)

rm(list=ls())

patients.matched = c("314", "1021247", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "A", "B", "C", "D", "X", "Y", "1009302", "1019357", "1020540", "1020583", "1021392", "1021631", "1022914", "1023056", "1023392", "1024518", "1024543", "1025108", "1025409", "1187", "FB11", "G", "HS6", "K", "MB2", "243", "933", "944", "KA14651", "KD20493", "KE12025", "MJ16441", "NH17331", "PJ13414", "RD17412", "RS13466", "ST13892", "ZA16211")

# get mutations
m <- data.frame(patient=character(), cohort=character(), sample=character(), gene=character(), chr=character(), pos=numeric(),	ref=character(), alt=character(), alt.reads=numeric(), tot.reads=numeric(), frequency=numeric())
reseq.relapsing <- read.delim("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", stringsAsFactors=F)

crebbp.dia <- reseq.relapsing[reseq.relapsing$status!="REJECT" & reseq.relapsing$non_silent==T & reseq.relapsing$sample=="rem_dia" & reseq.relapsing$gene=="CREBBP", c("patient", "sample", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
crebbp.dia$sample <- "diagnosis"
crebbp.dia <- data.frame(patient=crebbp.dia[,1], cohort="relapsing", crebbp.dia[,2:ncol(crebbp.dia)])
crebbp.dia$patient <- as.character(crebbp.dia$patient)
crebbp.dia$cohort <- as.character(crebbp.dia$cohort)
names(crebbp.dia) <- names(m)
m <- rbind(m, crebbp.dia)

# note: rel2 is actually first relapse for patient 715 (resequencing sample mix-up)
crebbp.rel <- reseq.relapsing[reseq.relapsing$status!="REJECT" & reseq.relapsing$non_silent==T & ((reseq.relapsing$sample=="rem_rel" & reseq.relapsing$patient != "715") | (reseq.relapsing$sample=="rem_rel2" & reseq.relapsing$patient == "715")) & reseq.relapsing$gene=="CREBBP", c("patient", "sample", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
crebbp.rel$sample <- "relapse"
crebbp.rel <- data.frame(patient=crebbp.rel[,1], cohort="relapsing", crebbp.rel[,2:ncol(crebbp.rel)])
crebbp.rel$patient <- as.character(crebbp.rel$patient)
crebbp.rel$cohort <- as.character(crebbp.rel$cohort)
names(crebbp.rel) <- names(m)
m <- rbind(m, crebbp.rel)

#cleanup
m$patient <- as.factor(m$patient)

# add columns
m$mut <- as.factor(paste0(m$gene, ":", m$chr, ":", m$pos, ":", m$ref, ">", m$alt))
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
m$frequency.norm[m$frequency.norm>1] <- m$frequency[m$frequency.norm>1]
m$patient.label <- as.factor(ifelse(m$patient %in% patients.matched, paste0(m$patient,"*"), as.character(m$patient)))

# manually add low AF conserved mutations identified by backtracking (= visual examination in GenomeBrowse)
r <- m[m$patient=="A",]; r$sample="diagnosis"; r$frequency.norm=0.02; m <- rbind(m,r)
r <- m[m$patient=="842",]; r$sample="diagnosis"; r$frequency.norm=0.02; m <- rbind(m,r)
r <- m[m$patient=="G",]; r$sample="diagnosis"; r$frequency.norm=0.02; m <- rbind(m,r)

# subset samples
m.dia <- m[m$sample=="diagnosis",]
m.dia <- ddply(m.dia[order(m.dia$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.rel <- m[m$sample=="relapse",]
m.rel <- ddply(m.rel[order(m.rel$frequency.norm, decreasing=T),], .(patient), transform, midpoint=cumsum(frequency.norm)-0.5*frequency.norm+0.005) # add midpoint for label plotting
m.cons <- merge(m.dia, m.rel, by=c("patient", "mut"))[,c("patient", "mut", "frequency.norm.x", "frequency.norm.y")]
m.cons <- paste0(m.cons$patient, ":", m.cons$mut)

pdf("~/hdall/results/crebbp-heterogeneity/crebbp-barplot.pdf", width=18, height=13)

#---
# RELAPSING, DIA
#---

# reorder patients by cumulative frequency
sorted <- aggregate(frequency.norm~patient.label, data=m.dia, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient.label", "frequency.norm")]
m.dia$patient.label <- factor(m.dia$patient.label, sorted$patient.label)

# find heterogenious samples
het.dia <- aggregate(frequency.norm~patient, data=m.dia, FUN=length)
het.dia.patient <- het.dia[het.dia$frequency.norm>1,"patient"]
#m.dia$group <- ifelse(m.dia$patient %in% het.dia.patient, sprintf("heterogenious (%.0f%%)", length(het.dia.patient)/length(levels(m.dia$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.dia.patient)/length(levels(m.dia$patient))*100))
m.dia$group <- ifelse(m.dia$patient %in% het.dia.patient, "heterogeneous", "homogeneous")
m.dia$label <- ifelse(paste0(m.dia$patient, ":", m.dia$mut) %in% m.cons, "symbol(\"\\267\")", NA)

plot.dia <- ggplot(data=m.dia, aes(x=patient.label, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(.~group, scale="free_x", space = "free_x") +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
#		geom_text(aes(label=sprintf("%.2f", frequency.norm), y=midpoint), size=4, vjust=1, colour="white", parse=TRUE) +
		scale_y_continuous(expand = c(0,0), limits = c(0,1.5), breaks=c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.position="none",
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
sorted <- aggregate(frequency.norm~patient.label, data=m.rel, FUN=sum)
sorted <- sorted[order(sorted$frequency.norm, decreasing=T), c("patient.label", "frequency.norm")]
m.rel$patient.label <- factor(m.rel$patient.label, sorted$patient.label)

# find heterogenious samples
het.rel <- aggregate(frequency.norm~patient, data=m.rel, FUN=length)
het.rel.patient <- het.rel[het.rel$frequency.norm>1,"patient"]
#m.rel$group <- ifelse(m.rel$patient %in% het.rel.patient, sprintf("heterogenious (%.0f%%)", length(het.rel.patient)/length(levels(m.rel$patient))*100), sprintf("homogenious (%.0f%%)", 100-length(het.rel.patient)/length(levels(m.rel$patient))*100))
m.rel$group <- ifelse(m.rel$patient %in% het.rel.patient, "heterogeneous", "homogeneous")
m.rel$label <- ifelse(paste0(m.rel$patient, ":", m.rel$mut) %in% m.cons, "symbol(\"\\267\")", NA)

plot.rel <- ggplot(data=m.rel, aes(x=patient.label, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
				facet_grid(.~group, scale="free_x", space = "free_x") +
				geom_bar(stat="identity", width=0.9, colour="white") +
				geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
#				geom_text(aes(label=sprintf("%.2f", frequency.norm), y=midpoint), size=4, vjust=1, colour="white", parse=TRUE) +
				scale_y_continuous(expand = c(0,0), limits = c(0,1.5), breaks=c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) +
				theme_bw() + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.x = element_text(size=20, vjust=0.2), axis.title.y = element_text(size=18, vjust=0.1), 
						legend.position="none",
						strip.text.x = element_text(size=15), 
						panel.grid.major.x = element_blank(),
						panel.grid.minor.y = element_blank(),
						plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(1,1,0.3,1), "cm")) +
#				geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
				ylab("Allelic frequency relapse") + 
				xlab("Case") 
		#				+ ggtitle("Ras pathway mutations at relapse")
		
grid.arrange(plot.dia, plot.rel, nrow=2)

#---
# P-value
#---
sig <- fisher.test(matrix(c(sum(het.dia$frequency.norm>1), sum(het.dia$frequency.norm==1), sum(het.rel$frequency.norm>1), sum(het.rel$frequency.norm==1)), nrow=2, byrow=T))
print(sprintf("P-value (Fisher's exact test): %g", sig$p.value))

dev.off()

#---
# RELAPSING, MATCHED, DIA+REL
#---

# reorder patients by cumulative frequency
m.combined <- m.dia[m.dia$patient %in% patients.matched,]
m.combined <- rbind(m.combined, m.rel[m.rel$patient %in% patients.matched,])

# reorder patients by group and cumulative frequency
sorted <- aggregate(frequency.norm~patient+sample+group, data=m.combined, FUN=sum)
sorted$sample <- as.character(sorted$sample)
sorted$sample[sorted$sample=="relapse"] <- "0relapse"
sorted$group[sorted$group=="homogeneous"] <- "0homogeneous"
sorted <- sorted[order(sorted$sample, sorted$group, sorted$frequency.norm, decreasing=T),]
m.combined$patient <- factor(m.combined$patient, unique(as.character(sorted$patient)))

pdf("~/hdall/results/crebbp-heterogeneity/crebbp-barplot-matched.pdf", width=14, height=9)
plot.dia <- ggplot(data=m.combined, aes(x=patient, y=frequency.norm, fill=mut, order=-frequency.norm)) + 
		facet_grid(sample~.) +
		geom_bar(stat="identity", width=0.9, colour="white") +
		geom_text(aes(label=label, y=midpoint), size=8, vjust=1, colour="white", parse=TRUE) +
#		geom_text(aes(label=sprintf("%.2f", frequency.norm), y=midpoint), size=4, vjust=1, colour="white", parse=TRUE) +
		scale_y_continuous(expand = c(0,0), limits = c(0,1.5), breaks=c(0.25, 0.5, 0.75, 1, 1.25, 1.5)) +
		theme_bw() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15), axis.text.y = element_text(size=13), axis.title.y = element_text(size=18, vjust=0.1), axis.title.x = element_blank(), 
				legend.position="none",
				strip.text.x = element_text(size=15), 
				panel.grid.major.x = element_blank(),
				panel.grid.minor.y = element_blank(),
				plot.title=element_text(size=20, face="bold", vjust=1), plot.margin=unit(c(0.5,1,0,1), "cm")) +
#		geom_text(aes(label = mut.short), size = 3, hjust = 0.5, vjust = 1.2, position = "stack", colour="white") +
		ylab("Allelic frequency")  
#		+ ggtitle("Ras pathway mutations at diagnosis")
print(plot.dia)
dev.off()

