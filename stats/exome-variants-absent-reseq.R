options(warn=1)

m <- read.delim("/mnt/projects/hdall/results/reseq/filtered-variants.merged.tsv", check.names=F, stringsAsFactor=F)
m <- m[m$sample != "rem_rel2" & m$patient != "E" & m$chr != "chr6_mcf_hap5",] 
g <- read.delim("/mnt/projects/hdall/results/panel-genes-paper.tsv", stringsAsFactor=F, header=F)

e <- m[!is.na(m$status.dis) & m$status.dis=="PASS" & m$non_silent.dis==1 & m$freq_leu.dis>=0.1 & m$gene %in% g$V1,]
r <- m[!is.na(m$status.val) & m$status.val=="PASS" & m$non_silent.val==1 & m$patient %in% e$patient,]

es <- e[(is.na(e$status.val) | e$status.val!="PASS") & e$patient != "715",] # exome-specific; exclude patient 715 because samples not matched b/w exome and reseq
rs <- r[(is.na(r$status.dis) | r$status.dis!="PASS" | r$freq_leu.dis < 0.1) & r$patient != "715",] # relapse-specific; exclude patient 715 because samples not matched b/w exome and reseq

print(sprintf("Number of non-silent exome sequencing variants with AF >= 10%% detected in panel genes: %d", nrow(e)))
print(sprintf("  not detected by resequencing: %d (%.1f%%)", nrow(es),nrow(es)/nrow(e)*100))
print("Summary of AFs of variants missed by resequencing")
print(summary(es$freq_leu.dis))

print("Missed by patient and sample")
es$patient <- as.factor(es$patient)
es$sample <- as.factor(es$sample)
print(table(es$patient, es$sample))

print(sprintf("Number of non-silent resequencing variants detected in exome sequencing patients: %d", nrow(r)))
print(sprintf("  not detected by exome sequencing: %d (%.1f%%)", nrow(rs),nrow(rs)/nrow(r)*100))

