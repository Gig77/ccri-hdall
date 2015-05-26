options(warn=1)

m <- read.delim("/mnt/projects/hdall/results/impacted-genes-list.tsv")
m$chr <- gsub("chr", "hs", as.character(m$chr))
m$loc <- paste(m$chr, m$start, m$end, m$gene, sep=":")

all <- aggregate(num_mut~loc, data=m[m$comparison=="rem_rel",], FUN=sum)
all <- read.table(text=paste(all$loc, all$num_mut, sep=":"), sep=":",col.names=c("chr", "start", "end", "gene", "num_mut"))
write.table(all[!is.na(all$chr),c("chr", "start", "end", "num_mut")], file="/mnt/projects/hdall/results/mutations-exome.all.circos.tsv", row.names=F, col.names=F, sep="\t", quote=F)

ns <- aggregate(num_mut_nonsyn~loc, data=m[m$comparison=="rem_rel",], FUN=sum)
ns <- read.table(text=paste(ns$loc, ns$num_mut_nonsyn, sep=":"), sep=":",col.names=c("chr", "start", "end", "gene", "num_mut"))
write.table(ns[!is.na(ns$num_mut), c("chr", "start", "end", "num_mut")], file="/mnt/projects/hdall/results/mutations-exome.nonsilent.circos.tsv", row.names=F, col.names=F, sep="\t", quote=F)

del <- aggregate(num_mut_deleterious~loc, data=m[m$comparison=="rem_rel",], FUN=sum)
del <- read.table(text=paste(del$loc, del$num_mut_deleterious, sep=":"), sep=":",col.names=c("chr", "start", "end", "gene", "num_mut"))
write.table(del[!is.na(del$num_mut),c("chr", "start", "end", "num_mut")], file="/mnt/projects/hdall/results/mutations-exome.deleterious.circos.tsv", row.names=F, col.names=F, sep="\t", quote=F)
		