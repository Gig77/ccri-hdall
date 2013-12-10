data <- read.csv("../filtered-variants.cosmic.normaf.tsv", sep="\t")
data <- data[data$status!="REJECT" & data$non_silent==1,]

dia <- data[data$sample=="rem_dia", c("patient", "chr", "pos", "var_type", "gene", "freq_leu_norm")]
rel <- data[data$sample=="rem_rel", c("patient", "chr", "pos", "var_type", "gene", "freq_leu_norm")]
names(dia)[6] <- "dia"
names(rel)[6] <- "rel"
m <- merge(dia, rel, all.x=T, all.y=T)

inc <- m[(is.na(m$dia) | m$dia<0.1) & (!is.na(m$rel) & m$rel>=0.2),]
dec <- m[(!is.na(m$dia) & m$dia>=0.2) & (is.na(m$rel) | m$rel<0.1),]
stable <- m[(!is.na(m$dia) & m$dia>=0.2) & (!is.na(m$rel) & m$rel>=0.2),]

inc.genes <- unique(as.character(inc$gene)) 
dec.genes <- unique(as.character(dec$gene)) 
stable.genes <- unique(as.character(stable$gene)) 

write(inc.genes, "genes.inc.txt", sep="\n")
write(dec.genes, "genes.dec.txt", sep="\n")
write(stable.genes, "genes.stable.txt", sep="\n")
