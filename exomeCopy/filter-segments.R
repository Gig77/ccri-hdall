s <- read.delim("~/hdall/results/exomeCopy/allpatients.compiled-segments.exomeCopy.tsv")
sf <- s
sf <- sf[sf$copy.count != 2,]
sf <- sf[!grepl("_rem", sf$overlap.samples),]

high.score <- sf$log.odds >= 50
conserved <- apply(sf, 1, function(x) { !is.na(x["overlap.samples"]) & grepl(sub("(_rem|_dia|_rel\\d?)", "", x["sample.name"], perl=T), x["overlap.samples"]) })
recurrent <- !is.na(sf$overlap.count.tumor) & sf$overlap.count.tumor >= 3

sf <- sf[high.score | (sf$nranges >= 2 & (conserved | recurrent)),]

write.table(sf, "~/hdall/results/exomeCopy/allpatients.filtered-segments.exomeCopy.tsv.part", row.names=F, col.names=T, sep="\t", quote=F)