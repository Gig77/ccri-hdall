rm(list=ls())

effect.cds <- c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT", "NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR", "START_GAINED", "START_LOST", "STOP_GAINED", "SYNONYMOUS_CODING", "SYNONYMOUS_STOP")
effect.nonsyn <- c("NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", "START_LOST", "STOP_GAINED")
effect.indel <- c("CODON_CHANGE_PLUS_CODON_DELETION", "CODON_CHANGE_PLUS_CODON_INSERTION", "CODON_DELETION", "CODON_INSERTION", "FRAME_SHIFT")
effect.splice <- c("SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR")
min.af <- 0.1

# exome seq variants
exome <- read.delim("/mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv", stringsAsFactor=F)
exome <- exome[exome$patient != "E",]
exome$patient <- as.factor(exome$patient)
reseq <- read.delim("/mnt/projects/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", stringsAsFactor=F)
reseq <- reseq[reseq$patient != "E",]
reseq$patient <- as.factor(reseq$patient)

exome.pass <- exome[exome$status != "REJECT",]
reseq.pass <- reseq[reseq$status != "REJECT",]

exome.pass.minaf <- exome.pass[exome.pass$freq_leu >= min.af,]
reseq.pass.minaf <- reseq.pass[reseq.pass$freq_leu >= min.af,]

exome.pass.minaf.dia <- exome.pass.minaf[exome.pass.minaf$sample=="rem_dia",]
exome.pass.minaf.rel <- exome.pass.minaf[exome.pass.minaf$sample=="rem_rel",]

merge.columns <- names(exome.pass.minaf)[1:17][!names(exome.pass.minaf)[1:17] %in% c("sample")] 
exome.pass.minaf.combined <- merge(exome.pass.minaf.dia, exome.pass.minaf.rel, by=merge.columns, all.x=T, all.y=T, suffixes=c(".dia", ".rel"))[,c(merge.columns, "freq_leu.dia", "freq_leu.rel")]

exome.pass.minaf.dia.cds <- exome.pass.minaf.dia[exome.pass.minaf.dia$effect %in% effect.cds,]
exome.pass.minaf.rel.cds <- exome.pass.minaf.rel[exome.pass.minaf.rel$effect %in% effect.cds,]
exome.pass.minaf.combined.cds <- exome.pass.minaf.combined[exome.pass.minaf.combined$effect %in% effect.cds,]

exome.pass.minaf.dia.cds.ns <- exome.pass.minaf.dia.cds[exome.pass.minaf.dia.cds$non_silent == 1,]
exome.pass.minaf.rel.cds.ns <- exome.pass.minaf.rel.cds[exome.pass.minaf.rel.cds$non_silent == 1,]
exome.pass.minaf.combined.cds.ns <- exome.pass.minaf.combined.cds[exome.pass.minaf.combined.cds$non_silent == 1,]

exome.pass.minaf.combined.cds.ns.nonsyn <- exome.pass.minaf.combined.cds.ns[exome.pass.minaf.combined.cds.ns$effect %in% effect.nonsyn,]
exome.pass.minaf.combined.cds.ns.indel <- exome.pass.minaf.combined.cds.ns[exome.pass.minaf.combined.cds.ns$effect %in% effect.indel,]
exome.pass.minaf.combined.cds.ns.splice <- exome.pass.minaf.combined.cds.ns[exome.pass.minaf.combined.cds.ns$effect %in% effect.splice,]

exome.pass.minaf.dia.cds.ns.del <- exome.pass.minaf.dia.cds.ns[exome.pass.minaf.dia.cds.ns$deleterious == "yes",]
exome.pass.minaf.rel.cds.ns.del <- exome.pass.minaf.rel.cds.ns[exome.pass.minaf.rel.cds.ns$deleterious == "yes",]
exome.pass.minaf.combined.cds.ns.del <- exome.pass.minaf.combined.cds.ns[exome.pass.minaf.combined.cds.ns$deleterious == "yes",]

exome.pass.minaf.dia.cds.ns.del.point <- exome.pass.minaf.dia.cds.ns.del[exome.pass.minaf.dia.cds.ns.del$var_type == "snp",]
exome.pass.minaf.rel.cds.ns.del.point <- exome.pass.minaf.rel.cds.ns.del[exome.pass.minaf.rel.cds.ns.del$var_type == "snp",]
exome.pass.minaf.combined.cds.ns.del.point <- exome.pass.minaf.combined.cds.ns.del[exome.pass.minaf.combined.cds.ns.del$var_type == "snp",]

exome.pass.minaf.dia.cds.ns.del.indel <- exome.pass.minaf.dia.cds.ns.del[exome.pass.minaf.dia.cds.ns.del$var_type == "indel",]
exome.pass.minaf.rel.cds.ns.del.indel <- exome.pass.minaf.rel.cds.ns.del[exome.pass.minaf.rel.cds.ns.del$var_type == "indel",]
exome.pass.minaf.combined.cds.ns.del.indel <- exome.pass.minaf.combined.cds.ns.del[exome.pass.minaf.combined.cds.ns.del$var_type == "indel",]

print(sprintf("exome -- AF>=%.1f -- diagnosis -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia), mean(table(factor(exome.pass.minaf.dia[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia[,"patient"])), max(table(exome.pass.minaf.dia[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- diagnosis -- CDS -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia.cds), mean(table(factor(exome.pass.minaf.dia.cds[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia.cds[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia.cds[,"patient"])), max(table(exome.pass.minaf.dia.cds[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- diagnosis -- CDS -- non-silent -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia.cds.ns), mean(table(factor(exome.pass.minaf.dia.cds.ns[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia.cds.ns[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia.cds.ns[,"patient"])), max(table(exome.pass.minaf.dia.cds.ns[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- diagnosis -- CDS -- non-silent -- deleterious -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia.cds.ns.del), mean(table(factor(exome.pass.minaf.dia.cds.ns.del[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia.cds.ns.del[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia.cds.ns.del[,"patient"])), max(table(exome.pass.minaf.dia.cds.ns.del[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- diagnosis -- CDS -- non-silent -- deleterious -- point mutations -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia.cds.ns.del.point), mean(table(factor(exome.pass.minaf.dia.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia.cds.ns.del.point[,"patient"])), max(table(exome.pass.minaf.dia.cds.ns.del.point[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- diagnosis -- CDS -- non-silent -- deleterious -- indels -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.dia.cds.ns.del.indel), mean(table(factor(exome.pass.minaf.dia.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.dia.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.dia.cds.ns.del.indel[,"patient"])), max(table(exome.pass.minaf.dia.cds.ns.del.indel[,"patient"]))))

print(sprintf("exome -- AF>=%.1f -- relapse -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel), mean(table(factor(exome.pass.minaf.rel[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel[,"patient"])), max(table(exome.pass.minaf.rel[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- relapse -- CDS -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel.cds), mean(table(factor(exome.pass.minaf.rel.cds[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel.cds[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel.cds[,"patient"])), max(table(exome.pass.minaf.rel.cds[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- relapse -- CDS -- non-silent -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel.cds.ns), mean(table(factor(exome.pass.minaf.rel.cds.ns[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel.cds.ns[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel.cds.ns[,"patient"])), max(table(exome.pass.minaf.rel.cds.ns[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- relapse -- CDS -- non-silent -- deleterious -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel.cds.ns.del), mean(table(factor(exome.pass.minaf.rel.cds.ns.del[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel.cds.ns.del[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel.cds.ns.del[,"patient"])), max(table(exome.pass.minaf.rel.cds.ns.del[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- relapse -- CDS -- non-silent -- deleterious -- point mutations -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel.cds.ns.del.point), mean(table(factor(exome.pass.minaf.rel.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel.cds.ns.del.point[,"patient"])), max(table(exome.pass.minaf.rel.cds.ns.del.point[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- relapse -- CDS -- non-silent -- deleterious -- indels -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.rel.cds.ns.del.indel), mean(table(factor(exome.pass.minaf.rel.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.rel.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.rel.cds.ns.del.indel[,"patient"])), max(table(exome.pass.minaf.rel.cds.ns.del.indel[,"patient"]))))

print(sprintf("exome -- AF>=%.1f -- combined -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined), mean(table(factor(exome.pass.minaf.combined[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined[,"patient"])), max(table(exome.pass.minaf.combined[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds), mean(table(factor(exome.pass.minaf.combined.cds[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds[,"patient"])), max(table(exome.pass.minaf.combined.cds[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns), mean(table(factor(exome.pass.minaf.combined.cds.ns[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- non-synonymous -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.nonsyn), mean(table(factor(exome.pass.minaf.combined.cds.ns.nonsyn[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.nonsyn[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.nonsyn[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.nonsyn[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- indels -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.indel), mean(table(factor(exome.pass.minaf.combined.cds.ns.indel[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.indel[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.indel[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.indel[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- splice site -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.splice), mean(table(factor(exome.pass.minaf.combined.cds.ns.splice[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.splice[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.splice[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.splice[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- deleterious -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.del), mean(table(factor(exome.pass.minaf.combined.cds.ns.del[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.del[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.del[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.del[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- deleterious -- point mutations -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.del.point), mean(table(factor(exome.pass.minaf.combined.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.del.point[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.del.point[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.del.point[,"patient"]))))
print(sprintf("exome -- AF>=%.1f -- combined -- CDS -- non-silent -- deleterious -- indels -- No. (mean, median, min-max): %d (%.2f, %d, %d-%d)", min.af, nrow(exome.pass.minaf.combined.cds.ns.del.indel), mean(table(factor(exome.pass.minaf.combined.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), median(table(factor(exome.pass.minaf.combined.cds.ns.del.indel[,"patient"], levels=levels(exome$patient)))), min(table(exome.pass.minaf.combined.cds.ns.del.indel[,"patient"])), max(table(exome.pass.minaf.combined.cds.ns.del.indel[,"patient"]))))

# mutation profile
mp <- read.delim("/mnt/projects/hdall/results/stats/mutation-profile.data.af10.tsv", stringsAsFactor=F)
mp[,"transition"] <- mp$mutation == "A>G" | mp$mutation == "C>T"
agg <- aggregate(cbind(freq.dia, freq.rel)~patient+transition, sum, data=mp)
tstv <- data.frame(patient=agg[agg$transition==T, "patient"], tstv.dia=agg[agg$transition==T,"freq.dia"]/agg[agg$transition==F,"freq.dia"], tstv.rel=agg[agg$transition==T,"freq.rel"]/agg[agg$transition==F,"freq.rel"])
print("Transition/transversion ratios")
print(tstv)
print(sprintf("Average Ts/Tv ratio at diagnosis: %.2f", mean(tstv$tstv.dia)))
print(sprintf("Average Ts/Tv ratio at relapse: %.2f", mean(tstv$tstv.rel)))