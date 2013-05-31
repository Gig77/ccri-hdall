jacc <- function(a, b)
{
	i <- intersect(a, b)
	u <- union(a, b)
	length(i) / length(u)
}

t <- read.csv("~/hdall/results/music/rel/sm_pathways.annotated.tsv", sep="\t", as.is=c(1, 9))

t2 <- t[t$Class=="KEGG_PATHWAY" | t$Class == "BBID" | t$Class == "BIOCARTA",]

pairs <- t(combn(rownames(t2),2))

m <- matrix(nrow=nrow(t2), ncol=nrow(t2))
rownames(m) <- t2$Pathway
colnames(m) <- t2$Pathway

for (i in 1:nrow(pairs))
{
	index1 <- as.integer(pairs[i,1])
	index2 <- as.integer(pairs[i,2])
	
	p1 <- t$Pathway[index1]
	p2 <- t$Pathway[index2]
	
	g1 <- strsplit(t$Genes[index1], ",")[[1]]
	g2 <- strsplit(t$Genes[index2], ",")[[1]]
	
	j <- jacc(g1, g2)
	
	m[p1,p2] <- j		
	m[p2,p1] <- j		
	m[p1,p1] <- 1		
	m[p2,p2] <- 1		
}

# write distance matrix
write.csv(m, file="~/hdall/results/music/rel/sm_pathways.distance-matrix.tsv")

# hierarchical clustering
d <- as.dist(1-m)
hc <- hclust(d, method="average")
pdf("~/hdall/results/music/rel/sm_pathways.hclust.pdf")
plot(hc, hang=-1, cex=0.2, xlab="Pathways", main="Pathways clustered by shared genes")
dev.off()

# cut tree and output clusters
c <- cutree(hc, h=c(0.2,0.4,0.6,0.8))
c <- cbind(rownames(c), c)
colnames(c)[1] <- "Pathway"
merged <- merge(t, c, all.x=T)

write.csv(merged, file="~/hdall/results/music/rel/sm_pathways.annotated.clustered.tsv")


