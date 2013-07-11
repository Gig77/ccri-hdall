d <- read.csv("~/hdall/results/music/dia/sm_pathways.annotated.tsv", sep="\t", as.is=c(1, 9))
r <- read.csv("~/hdall/results/music/rel/sm_pathways.annotated.tsv", sep="\t", as.is=c(1, 9))
m <- merge(d, r, by=c("Pathway", "Name", "Class"), all=T, suffixes=c(".dia", ".rel"))
m <- m[order(m$p.value.rel),]
write.table(m, row.names=F, sep="\t")

