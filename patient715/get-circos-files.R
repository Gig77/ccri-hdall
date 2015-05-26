warnings()

library(vcd)
library(ggplot2)

min.af <- 0.10
cn.color <- data.frame(cnumber=c(0, 1, 2, 3, 4), color=c("fill_color=black", "fill_color=lred", "fill_color=lblue", "fill_color=lblue", "fill_color=dblue"))  # copy number 2 events are gains on chr X

# diagnosis, relapse 1
kamilla <- read.delim("/mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv", stringsAsFactors=F, na.strings=c("NA", "n/d"))
maria <- read.delim("/mnt/projects/p2ry8-crlf2/results/filtered-variants.cosmic.normaf.tsv", stringsAsFactor=F, na.strings=c("NA", "n/d"))
cnv <- read.delim("/mnt/projects/hdall/results/cnv/hdall.cnv.tsv", stringsAsFactor=F)
cnv <- cnv[cnv$patient=="715",]

m <- kamilla[kamilla$patient=="715" & kamilla$status != "REJECT" & kamilla$freq_leu >= min.af,]
m <- rbind(m, maria[maria$patient=="715" & maria$status != "REJECT" & maria$freq_leu >= min.af,])

m$chr <- gsub("chr", "hs", m$chr)
cnv$chromosome <- paste0("hs", cnv$chromosome)
cnv <- merge(cnv, cn.color, all.x=T)

m.dia <- m[m$sample=="rem_dia",]
m.rel <- m[m$sample=="rem_rel",]
m.rel3 <- m[m$sample=="rem_rel3",]

write.table(data.frame(m.dia[,"chr"], m.dia[,"pos"], m.dia[,"pos"]+1), file="/mnt/projects/hdall/results/patient715/mutations.dia.circos", col.names=F, row.names=F, sep="\t", quote=F)
write.table(data.frame(m.rel[,"chr"], m.rel[,"pos"], m.rel[,"pos"]+1), file="/mnt/projects/hdall/results/patient715/mutations.rel.circos", col.names=F, row.names=F, sep="\t", quote=F)
write.table(data.frame(m.rel3[,"chr"], m.rel3[,"pos"], m.rel3[,"pos"]+1), file="/mnt/projects/hdall/results/patient715/mutations.rel3.circos", col.names=F, row.names=F, sep="\t", quote=F)
write.table(cnv[cnv$sample=="rem_dia", c("chromosome", "start", "end", "color")], file="/mnt/projects/hdall/results/patient715/cnv.dia.circos", col.names=F, row.names=F, sep="\t", quote=F)
write.table(cnv[cnv$sample=="rem_rel", c("chromosome", "start", "end", "color")], file="/mnt/projects/hdall/results/patient715/cnv.rel.circos", col.names=F, row.names=F, sep="\t", quote=F)
write.table(cnv[cnv$sample=="rem_rel3", c("chromosome", "start", "end", "color")], file="/mnt/projects/hdall/results/patient715/cnv.rel3.circos", col.names=F, row.names=F, sep="\t", quote=F)
