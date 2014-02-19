library("RColorBrewer")
warnings()

rm(list=ls())

patients <- c("1021247", "446", "818", "A", "592", "430", "792", "786", "B", "460", "545", "715", "X", "Y", "842", "C", "D", "399", "314") 
#patients <- c("314", "399") 

# TABLE: filtered-variants.tsv
t <- read.delim("filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$var_type == "snp",]
min.af <- 0.10
max.af <- 1.00
min.dp.leu <- 10
min.dp.rem <- 10
color <- "Set1"

df <- data.frame(patient=character(), mutation=character(), freq.dia=numeric(), freq.rel=numeric(), stringsAsFactors=F)
m.dia <- matrix(rep(0, 6*length(patients)), nrow=length(patients))
colnames(m.dia) <- c("A>G", "C>T", "G>T", "A>C", "G>C", "A>T")
rownames(m.dia) <- patients
m.rel <- m.dia

# per patient

#pdf("stats/mutation-profile-stacked.pdf", width=15)
png("stats/mutation-profile-stacked.png")

#layout(matrix(c(seq(1:20),rep(21,5)), ncol=5, byrow=T), heights=c(rep(0.23, 4), 0.08))
#par(mar=c(2.0, 2.5, 3, 0))
for(p in patients) {
	# diagnosis
	tdia <- t[t$patient==p & t$sample=="rem_dia" & t$freq_leu>=min.af & t$freq_leu<=max.af & t$dp_rem_tot>=min.dp.rem & t$dp_leu_tot>=min.dp.leu,]
	m.dia[p,1] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="G")|(tdia$ref=="G" & tdia$alt=="A"),])
	m.dia[p,2] <- nrow(tdia[(tdia$ref=="C" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="C"),])
	m.dia[p,3] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="G"),])
	m.dia[p,4] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="C")|(tdia$ref=="C" & tdia$alt=="A"),])
	m.dia[p,5] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="C")|(tdia$ref=="C" & tdia$alt=="G"),])
	m.dia[p,6] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="A"),])
	m.dia[p,] <- m.dia[p,]/sum(m.dia[p,])
	
	# relapse
	trel <- t[t$patient==p & t$sample=="rem_rel" & t$freq_leu>=min.af & t$freq_leu<=max.af & t$dp_rem_tot>=min.dp.rem & t$dp_leu_tot>=min.dp.leu,]
	m.rel[p,1] <- nrow(trel[(trel$ref=="A" & trel$alt=="G")|(trel$ref=="G" & trel$alt=="A"),])
	m.rel[p,2] <- nrow(trel[(trel$ref=="C" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="C"),])
	m.rel[p,3] <- nrow(trel[(trel$ref=="G" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="G"),])
	m.rel[p,4] <- nrow(trel[(trel$ref=="A" & trel$alt=="C")|(trel$ref=="C" & trel$alt=="A"),])
	m.rel[p,5] <- nrow(trel[(trel$ref=="G" & trel$alt=="C")|(trel$ref=="C" & trel$alt=="G"),])
	m.rel[p,6] <- nrow(trel[(trel$ref=="A" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="A"),])
	m.rel[p,] <- m.rel[p,]/sum(m.rel[p,])

	df[nrow(df)+1,] <- c(p, "A>G", m.dia[p,1], m.rel[p,1])
	df[nrow(df)+1,] <- c(p, "C>T", m.dia[p,2], m.rel[p,2])
	df[nrow(df)+1,] <- c(p, "G>T", m.dia[p,3], m.rel[p,3])
	df[nrow(df)+1,] <- c(p, "A>C", m.dia[p,4], m.rel[p,4])
	df[nrow(df)+1,] <- c(p, "G>C", m.dia[p,5], m.rel[p,5])
	df[nrow(df)+1,] <- c(p, "A>T", m.dia[p,6], m.rel[p,6])
}

layout(matrix(c(1,2,3,3), nrow = 2), widths = c(0.8, 0.2))
par()
barplot(t(m.dia), beside=FALSE, col=brewer.pal(6, color), las=2, main="diagnosis")
par()
barplot(t(m.rel), beside=FALSE, col=brewer.pal(6, color), las=2, main="relapse")
par(mar=c(0, 0, 4, 0))
plot(0:1, 0:1, type="n", axes=F, ann=F)
legend("topleft", rev(colnames(m.dia)), fill=rev(brewer.pal(6, color)))

write.table(df, file="stats/mutation-profile.data.af10.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#par(xpd=TRUE, mai=c(0,0,0,0))
#plot.new()
#legend(x="center", legend=c("diagnosis", "relapse"), fill=c("black", "white"), horiz=T)

dev.off()
