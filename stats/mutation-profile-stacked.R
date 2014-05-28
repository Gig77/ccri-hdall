library("RColorBrewer")
library(vcd)
warnings()

rm(list=ls())

patients <- c("1021247", "446", "818", "A", "592", "430", "792", "786", "B", "460", "545", "715", "X", "Y", "842", "C", "D", "399", "314") 
#patients <- c("314", "399") 
min.af <- 0.10
max.af <- 1.00
min.dp.leu <- 10
min.dp.rem <- 10
color <- "Set1"

# TABLE: filtered-variants.tsv
t <- read.delim("~/hdall/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$var_type == "snp",]

t$class <- NA
t$class[(t$ref=="A" & t$alt=="G")|(t$ref=="T" & t$alt=="C")] <- "A:T>G:C" # transition
t$class[(t$ref=="G" & t$alt=="T")|(t$ref=="C" & t$alt=="A")] <- "G:C>T:A" # transvertion
t$class[(t$ref=="G" & t$alt=="A")|(t$ref=="C" & t$alt=="T")] <- "G:C>A:T" # transition
t$class[(t$ref=="A" & t$alt=="C")|(t$ref=="T" & t$alt=="G")] <- "A:T>C:G" # transvertion
t$class[(t$ref=="G" & t$alt=="C")|(t$ref=="C" & t$alt=="G")] <- "G:C>C:G" # transvertion
t$class[(t$ref=="A" & t$alt=="T")|(t$ref=="T" & t$alt=="A")] <- "A:T>T:A" # transvertion

t.dia <- t[t$sample=="rem_dia", c("patient", "chr", "pos", "ref", "alt", "gene", "class", "freq_leu", "dp_rem_tot", "dp_leu_tot")]
names(t.dia)[8:10] <- c("freq_leu.dia", "dp_rem_tot.dia", "dp_leu_tot.dia")
t.rel <- t[t$sample=="rem_rel", c("patient", "chr", "pos", "ref", "alt", "gene", "class", "freq_leu", "dp_rem_tot", "dp_leu_tot")]
names(t.rel)[8:10] <- c("freq_leu.rel", "dp_rem_tot.rel", "dp_leu_tot.rel")
t.rel.specific <- merge(t.dia, t.rel, all.x=T, all.y=T)

t.dia <- t.dia[t.dia$freq_leu.dia>=min.af & t.dia$freq_leu.dia<=max.af & t.dia$dp_rem_tot.dia>=min.dp.rem & t.dia$dp_leu_tot.dia>=min.dp.leu,]
t.rel <- t.rel[t.rel$freq_leu.rel>=min.af & t.rel$freq_leu.rel<=max.af & t.rel$dp_rem_tot.rel>=min.dp.rem & t.rel$dp_leu_tot.rel>=min.dp.leu,]
t.rel.specific <- t.rel.specific[is.na(t.rel.specific$freq_leu.dia) & t.rel.specific$freq_leu.rel>=min.af & t.rel.specific$freq_leu.rel<=max.af & t.rel.specific$dp_rem_tot.rel>=min.dp.rem & t.rel.specific$dp_leu_tot.rel>=min.dp.leu,]

df <- data.frame(patient=character(), mutation=character(), freq.dia=numeric(), freq.rel=numeric(), freq.rel.specific=numeric(), stringsAsFactors=F)
m.dia <- matrix(rep(0, 6*(length(patients)+1)), nrow=length(patients)+1)
colnames(m.dia) <- c("G:C>A:T", "G:C>T:A", "A:T>G:C", "A:T>C:G", "G:C>C:G", "A:T>T:A")
rownames(m.dia) <- c(patients, "overall")
m.rel <- m.dia
m.rel.specific <- m.dia

# per patient
#-------------

#layout(matrix(c(seq(1:20),rep(21,5)), ncol=5, byrow=T), heights=c(rep(0.23, 4), 0.08))
#par(mar=c(2.0, 2.5, 3, 0))
for(p in patients) {

	# diagnosis
	tdia <- t.dia[t.dia$patient==p,]
	m.dia[p,1] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="A")|(tdia$ref=="C" & tdia$alt=="T"),])
	m.dia[p,2] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="T")|(tdia$ref=="C" & tdia$alt=="A"),])
	m.dia[p,3] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="G")|(tdia$ref=="T" & tdia$alt=="C"),])
	m.dia[p,4] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="C")|(tdia$ref=="T" & tdia$alt=="G"),])
	m.dia[p,5] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="C")|(tdia$ref=="C" & tdia$alt=="G"),])
	m.dia[p,6] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="A"),])
	m.dia[p,] <- m.dia[p,]/sum(m.dia[p,])
	
	# relapse
	trel <- t.rel[t.rel$patient==p,]
	m.rel[p,1] <- nrow(trel[(trel$ref=="G" & trel$alt=="A")|(trel$ref=="C" & trel$alt=="T"),])
	m.rel[p,2] <- nrow(trel[(trel$ref=="G" & trel$alt=="T")|(trel$ref=="C" & trel$alt=="A"),])
	m.rel[p,3] <- nrow(trel[(trel$ref=="A" & trel$alt=="G")|(trel$ref=="T" & trel$alt=="C"),])
	m.rel[p,4] <- nrow(trel[(trel$ref=="A" & trel$alt=="C")|(trel$ref=="T" & trel$alt=="G"),])
	m.rel[p,5] <- nrow(trel[(trel$ref=="G" & trel$alt=="C")|(trel$ref=="C" & trel$alt=="G"),])
	m.rel[p,6] <- nrow(trel[(trel$ref=="A" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="A"),])
	m.rel[p,] <- m.rel[p,]/sum(m.rel[p,])

	# relapse-specific
	trelspec <- t.rel.specific[t.rel.specific$patient==p,]
	m.rel.specific[p,1] <- nrow(trelspec[(trelspec$ref=="G" & trelspec$alt=="A")|(trelspec$ref=="C" & trelspec$alt=="T"),])
	m.rel.specific[p,2] <- nrow(trelspec[(trelspec$ref=="G" & trelspec$alt=="T")|(trelspec$ref=="C" & trelspec$alt=="A"),])
	m.rel.specific[p,3] <- nrow(trelspec[(trelspec$ref=="A" & trelspec$alt=="G")|(trelspec$ref=="T" & trelspec$alt=="C"),])
	m.rel.specific[p,4] <- nrow(trelspec[(trelspec$ref=="A" & trelspec$alt=="C")|(trelspec$ref=="T" & trelspec$alt=="G"),])
	m.rel.specific[p,5] <- nrow(trelspec[(trelspec$ref=="G" & trelspec$alt=="C")|(trelspec$ref=="C" & trelspec$alt=="G"),])
	m.rel.specific[p,6] <- nrow(trelspec[(trelspec$ref=="A" & trelspec$alt=="T")|(trelspec$ref=="T" & trelspec$alt=="A"),])
	m.rel.specific[p,] <- m.rel.specific[p,]/sum(m.rel.specific[p,])
	
	df[nrow(df)+1,] <- c(p, "G:C>A:T", m.dia[p,1], m.rel[p,1], m.rel.specific[p,1])
	df[nrow(df)+1,] <- c(p, "G:C>T:A", m.dia[p,2], m.rel[p,2], m.rel.specific[p,2])
	df[nrow(df)+1,] <- c(p, "A:T>G:C", m.dia[p,3], m.rel[p,3], m.rel.specific[p,3])
	df[nrow(df)+1,] <- c(p, "A:T>C:G", m.dia[p,4], m.rel[p,4], m.rel.specific[p,4])
	df[nrow(df)+1,] <- c(p, "G:C>C:G", m.dia[p,5], m.rel[p,5], m.rel.specific[p,5])
	df[nrow(df)+1,] <- c(p, "A:T>T:A", m.dia[p,6], m.rel[p,6], m.rel.specific[p,6])
}

# all patients combined
#-----------------------

# diagnosis
m.dia["overall",1] <- nrow(t.dia[(t.dia$ref=="G" & t.dia$alt=="A")|(t.dia$ref=="C" & t.dia$alt=="T"),])
m.dia["overall",2] <- nrow(t.dia[(t.dia$ref=="G" & t.dia$alt=="T")|(t.dia$ref=="C" & t.dia$alt=="A"),])
m.dia["overall",3] <- nrow(t.dia[(t.dia$ref=="A" & t.dia$alt=="G")|(t.dia$ref=="T" & t.dia$alt=="C"),])
m.dia["overall",4] <- nrow(t.dia[(t.dia$ref=="A" & t.dia$alt=="C")|(t.dia$ref=="T" & t.dia$alt=="G"),])
m.dia["overall",5] <- nrow(t.dia[(t.dia$ref=="G" & t.dia$alt=="C")|(t.dia$ref=="C" & t.dia$alt=="G"),])
m.dia["overall",6] <- nrow(t.dia[(t.dia$ref=="A" & t.dia$alt=="T")|(t.dia$ref=="T" & t.dia$alt=="A"),])
m.dia["overall",] <- m.dia["overall",]/sum(m.dia["overall",])

# relapse
m.rel["overall",1] <- nrow(t.rel[(t.rel$ref=="G" & t.rel$alt=="A")|(t.rel$ref=="C" & t.rel$alt=="T"),])
m.rel["overall",2] <- nrow(t.rel[(t.rel$ref=="G" & t.rel$alt=="T")|(t.rel$ref=="C" & t.rel$alt=="A"),])
m.rel["overall",3] <- nrow(t.rel[(t.rel$ref=="A" & t.rel$alt=="G")|(t.rel$ref=="T" & t.rel$alt=="C"),])
m.rel["overall",4] <- nrow(t.rel[(t.rel$ref=="A" & t.rel$alt=="C")|(t.rel$ref=="T" & t.rel$alt=="G"),])
m.rel["overall",5] <- nrow(t.rel[(t.rel$ref=="G" & t.rel$alt=="C")|(t.rel$ref=="C" & t.rel$alt=="G"),])
m.rel["overall",6] <- nrow(t.rel[(t.rel$ref=="A" & t.rel$alt=="T")|(t.rel$ref=="T" & t.rel$alt=="A"),])
m.rel["overall",] <- m.rel["overall",]/sum(m.rel["overall",])

# relapse-specific
m.rel.specific["overall",1] <- nrow(t.rel.specific[(t.rel.specific$ref=="G" & t.rel.specific$alt=="A")|(t.rel.specific$ref=="C" & t.rel.specific$alt=="T"),])
m.rel.specific["overall",2] <- nrow(t.rel.specific[(t.rel.specific$ref=="G" & t.rel.specific$alt=="T")|(t.rel.specific$ref=="C" & t.rel.specific$alt=="A"),])
m.rel.specific["overall",3] <- nrow(t.rel.specific[(t.rel.specific$ref=="A" & t.rel.specific$alt=="G")|(t.rel.specific$ref=="T" & t.rel.specific$alt=="C"),])
m.rel.specific["overall",4] <- nrow(t.rel.specific[(t.rel.specific$ref=="A" & t.rel.specific$alt=="C")|(t.rel.specific$ref=="T" & t.rel.specific$alt=="G"),])
m.rel.specific["overall",5] <- nrow(t.rel.specific[(t.rel.specific$ref=="G" & t.rel.specific$alt=="C")|(t.rel.specific$ref=="C" & t.rel.specific$alt=="G"),])
m.rel.specific["overall",6] <- nrow(t.rel.specific[(t.rel.specific$ref=="A" & t.rel.specific$alt=="T")|(t.rel.specific$ref=="T" & t.rel.specific$alt=="A"),])
m.rel.specific["overall",] <- m.rel.specific["overall",]/sum(m.rel.specific["overall",])


n.dia <- table(t.dia$patient)
n.dia[length(n.dia)+1] = sum(n.dia)
names(n.dia)[length(n.dia)]<-"overall"
for(i in 1:length(rownames(m.dia))) { rownames(m.dia)[i] <- paste0(rownames(m.dia)[i], " (n=", n.dia[rownames(m.dia)[i]], ")") }

n.rel <- table(t.rel$patient)
n.rel[length(n.rel)+1] = sum(n.rel)
names(n.rel)[length(n.rel)]<-"overall"
for(i in 1:length(rownames(m.rel))) { rownames(m.rel)[i] <- paste0(rownames(m.rel)[i], " (n=", n.rel[rownames(m.rel)[i]], ")") }

n.rel.specific <- table(t.rel.specific$patient)
n.rel.specific[length(n.rel.specific)+1] = sum(n.rel.specific)
names(n.rel.specific)[length(n.rel.specific)]<-"overall"
for(i in 1:length(rownames(m.rel.specific))) { rownames(m.rel.specific)[i] <- paste0(rownames(m.rel.specific)[i], " (n=", n.rel.specific[rownames(m.rel.specific)[i]], ")") }

do_plot <- function() {
#	layout(matrix(c(1,2,3,4,4,4), nrow = 3), widths = c(0.8, 0.2))
	layout(matrix(c(1,2,3,3), nrow = 2), widths = c(0.8, 0.2))
	par(mar=c(7,4,3,1))
	barplot(t(m.dia), beside=FALSE, col=brewer.pal(6, color), las=2, main="mutations at diagnosis", ylab="Proportion of substitutions (%)")
#	par()
#	barplot(t(m.rel), beside=FALSE, col=brewer.pal(6, color), las=2, main="relapse")
	par(mar=c(7,4,3,1))
	barplot(t(m.rel.specific), beside=FALSE, col=brewer.pal(6, color), las=2, main="relapse-specific mutations", ylim=c(0,1.035), ylab="Proportion of substitutions (%)")

	# indicate with asterix patients with significant difference in mutation profile b/w diagnosis and relapse
	text(2.0, 1.025, "*", cex=1.5) # 446
	text(5.55, 1.025, "*", cex=1.5) # 592
	text(11.6, 1.025, "*", cex=1.5) # 460
#	text(16.35, 1.025, "*", cex=1.5) # Y
	text(18.7, 1.025, "*", cex=1.5) # C
	text(19.9, 1.025, "*", cex=1.5) # D
#	text(21.1, 1.025, "*", cex=1.5) # 399
	text(23.5, 1.025, "**", cex=1.5) # all
	par(mar=c(0, 0, 4, 0))
	plot(0:1, 0:1, type="n", axes=F, ann=F)
	legend("topleft", rev(colnames(m.dia)), fill=rev(brewer.pal(6, color)))
}

#pdf("~/hdall/results/stats/mutation-profile-stacked.pdf", width=8, height=10)
pdf("~/hdall/results/stats/mutation-profile-stacked.pdf", width=8, height=8)
do_plot()
dev.off()

#png("~/hdall/results/stats/mutation-profile-stacked.png", width=2000, height=3000, res=200)
png("~/hdall/results/stats/mutation-profile-stacked.png", width=2000, height=2000, res=200)
do_plot()
dev.off()

print(sprintf("Average ratio Ts/Tv at diagnosis: %.3f", mean(aggregate(as.numeric(freq.dia) ~ patient, data = df[df$mutation %in% c("A:T>G:C", "G:C>A:T"),], sum)[,2]/aggregate(as.numeric(freq.dia) ~ patient, data = df[!df$mutation %in% c("A:T>G:C", "G:C>A:T"),], sum)[,2])))
print(sprintf("Average ratio Ts/Tv at relapse: %.3f", mean(aggregate(as.numeric(freq.rel.specific) ~ patient, data = df[df$mutation %in% c("A:T>G:C", "G:C>A:T"),], sum)[,2]/aggregate(as.numeric(freq.rel.specific) ~ patient, data = df[!df$mutation %in% c("A:T>G:C", "G:C>A:T"),], sum)[,2])))

write.table(df, file="~/hdall/results/stats/mutation-profile.data.af10.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# run statistical tests on different mutation profiles in diagnosis vs. relapse
t.filtered <- t.dia
names(t.filtered)[8:10] <- c("freq_leu", "dp_rem_tot", "dp_leu_tot")
t.filtered$sample <- "diagnosis"
tmp <- t.rel.specific[, c("patient", "chr", "pos", "ref", "alt", "gene", "class", "freq_leu.rel", "dp_rem_tot.rel", "dp_leu_tot.rel")]
names(tmp)[8:10] <- c("freq_leu", "dp_rem_tot", "dp_leu_tot")
tmp$sample <- "relapse"
t.filtered <- rbind(t.filtered, tmp)

pdf("~/hdall/results/stats/mutation-profile-stats-dia-vs-rel.pdf")
set.seed(22)
p <- fisher.test(table(t.filtered$sample, t.filtered$class), simulate.p.value=T,B=1000000)$p.value
tab <- mosaic(class ~ sample, data=t.filtered, pop=F, main=sprintf("mutations from all patients (n=%d, p=%.2g)", nrow(t.filtered), p))
labeling_cells(text=as.table(tab))(as.table(tab))
for(pat in patients)
{
	tfp <- t.filtered[t.filtered$patient==pat,]
	p <- fisher.test(table(tfp$sample, tfp$class), simulate.p.value=T,B=1000000)$p.value
	tab <- mosaic(class ~ sample, data=tfp, pop=F, main=sprintf("patient %s (n=%d, p=%.2g)", pat, nrow(tfp), p))
	labeling_cells(text=as.table(tab))(as.table(tab))
}
dev.off()

#par(xpd=TRUE, mai=c(0,0,0,0))
#plot.new()
#legend(x="center", legend=c("diagnosis", "relapse"), fill=c("black", "white"), horiz=T)
