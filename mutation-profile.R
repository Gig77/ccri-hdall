patients <- c("314", "399", "430", "446", "460", "545", "592", "715", "786", "792", "818", "842", "1021247", "A", "B", "C", "D", "E", "X", "Y") 
#patients <- c("314", "399") 

t <- read.delim("~/hdall/results/filtered-variants.tsv", stringsAsFactors=F)

m <- matrix(rep(0, 12), nrow=2)
colnames(m) <- c("A>T", "A>G", "A>C", "G>C", "G>T", "C>T")
rownames(m) <- c("diagnosis", "relapse")

pdf("~/hdall/results/stats/mutation-profile.pdf", width=15)

layout(matrix(c(seq(1:20),rep(21,5)), ncol=5, byrow=T), heights=c(rep(0.23, 4), 0.08))
par(mar=c(2.0, 2.5, 3, 0))
for(p in patients) {
	# diagnosis
	tdia <- t[t$patient==p & t$sample=="rem_dia",]
	m[1,1] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="A"),])
	m[1,2] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="G")|(tdia$ref=="G" & tdia$alt=="A"),])
	m[1,3] <- nrow(tdia[(tdia$ref=="A" & tdia$alt=="C")|(tdia$ref=="C" & tdia$alt=="A"),])
	m[1,4] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="C")|(tdia$ref=="C" & tdia$alt=="G"),])
	m[1,5] <- nrow(tdia[(tdia$ref=="G" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="G"),])
	m[1,6] <- nrow(tdia[(tdia$lref=="C" & tdia$alt=="T")|(tdia$ref=="T" & tdia$alt=="C"),])

	# relapse
	trel <- t[t$patient==p & t$sample=="rem_rel",]
	m[2,1] <- nrow(trel[(trel$ref=="A" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="A"),])
	m[2,2] <- nrow(trel[(trel$ref=="A" & trel$alt=="G")|(trel$ref=="G" & trel$alt=="A"),])
	m[2,3] <- nrow(trel[(trel$ref=="A" & trel$alt=="C")|(trel$ref=="C" & trel$alt=="A"),])
	m[2,4] <- nrow(trel[(trel$ref=="G" & trel$alt=="C")|(trel$ref=="C" & trel$alt=="G"),])
	m[2,5] <- nrow(trel[(trel$ref=="G" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="G"),])
	m[2,6] <- nrow(trel[(trel$ref=="C" & trel$alt=="T")|(trel$ref=="T" & trel$alt=="C"),])

	# compute relative frequency
	mrel <- m
	mrel[1,] = m[1,]/sum(m[1,])
	mrel[2,] = m[2,]/sum(m[2,])
	
	barplot(mrel, beside=T, legend=F, main=paste("patient ", p, " (n=", sum(m[1,]), "+", sum(m[2,]), ")", sep=""), ylim=c(0, 0.4), col=c("black", "white"))
}

par(xpd=TRUE, mai=c(0,0,0,0))
plot.new()
legend(x="center", legend=c("diagnosis", "relapse"), fill=c("black", "white"), horiz=T)

dev.off()