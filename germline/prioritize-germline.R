m <- read.delim("/mnt/projects/hdall/results/germline/filtered-variants.germline.nonsilent.tsv", na.strings=c(""), stringsAsFactors=F)
m <- m[m$patient != "E.rel.exome",] # exclude patient E
m <- m[m$type != "somatic",] # exclude somatic
m <- m[m$both_strands_rem == "yes" | m$both_strands_leu == "yes",] # require support from both strands
m$selected_because <- NA
m$selected_because[!is.na(m$clinvar_variant)] <- "ClinVar pathogenic"
m$selected_because[m$gene=="TP53"] <- "TP53"
m$selected_because[m$gene=="CREBBP"] <- "CREBBP"
m$selected_because[m$gene=="CREBBP" & m$deleterious=="yes"] <- "CREBBP, deleterious"
m$selected_because[m$gene=="CREBBP" & (m$cosmic_hits_nt > 0 | m$cosmic_hits_aa > 0)] <- "CREBBP, in COSMIC"
m$selected_because[is.na(m$selected_because) & m$type=="LOH" & m$deleterious=="yes" & m$pval_leu < 0.05] <- "LOH, deleterious"
m$selected_because[is.na(m$selected_because) & m$cosmic_hits_nt>=5] <- "reported >= 5 times in COSMIC"
m$selected_because[is.na(m$selected_because) & m$dbSNP=="." & m$cosmic_hits_nt > 0] <- "unknown or rare SNP, reported in COSMIC"
write.table(m[!is.na(m$selected_because),], file="/mnt/projects/hdall/results/germline/filtered-variants.germline.nonsilent.high-priority.tsv", col.names=T, row.names=F, sep="\t", quote=F)