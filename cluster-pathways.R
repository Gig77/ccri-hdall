options(warn=1)

jacc <- function(a, b)
{
	i <- intersect(a, b)
	u <- union(a, b)
	length(i) / length(u)
}

args <- commandArgs(trailingOnly = TRUE)
if (is.na(args[1])) stop("input file not specified")
if (is.na(args[2])) stop("output file not specified")

t <- read.csv(args[1], sep="\t", as.is=TRUE, check.names=F)
#t <- read.csv("~/hdall/results/music/sm_pathways.tsv.part", sep="\t", as.is=TRUE)
#t2 <- t

pcol <- as.numeric(args[3]) # number of column containing p-value for filtering
pval <- as.numeric(args[4]) # maximum p-value
if (!is.na(pcol))
{
	t2 <- t[(!is.na(t[pcol]) & t[pcol] <= pval) | t$Class == "NCI",]
} else
{
	t2 <- t
}
#t2 <- t[t$Class=="KEGG_PATHWAY",]

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
	
#	g1 <- strsplit(t$Genes.rel[index1], ",")[[1]]
#	g2 <- strsplit(t$Genes.rel[index2], ",")[[1]]
	g1 <- union(strsplit(t$Genes.dia[index1], ",")[[1]], strsplit(t$Genes.rel[index1], ",")[[1]])
	g2 <- union(strsplit(t$Genes.dia[index2], ",")[[1]], strsplit(t$Genes.rel[index2], ",")[[1]])
	
	j <- jacc(g1[!is.na(g1)], g2[!is.na(g2)])
#	j <- jacc(g1, g2)
	
	m[p1,p2] <- j		
	m[p2,p1] <- j		
	m[p1,p1] <- 1		
	m[p2,p2] <- 1		
}

# write distance matrix
#write.csv(m, file="~/hdall/results/music/rel/sm_pathways.distance-matrix.tsv")

# hierarchical clustering
d <- as.dist(1-m)
hc <- hclust(d, method="average")
#pdf("~/hdall/results/music/rel/sm_pathways.hclust.pdf", width=20)
#plot(hc, hang=-1, cex=0.05, xlab="Pathways", main="Pathways clustered by shared genes")
#dev.off()

# cut tree and output concatenated string of cluster ids for sorting
c <- cutree(hc, h=seq(0.18,0.98,by=0.2))
c <- c[hc$order,] # reorder pathways based on hierarchical tree
c <- cbind(1:nrow(c), c)
m <- as.data.frame(apply(matrix(sprintf("%04d", c), ncol=6), 1, paste, collapse="-"))
m <- cbind(rownames(c), m)
colnames(m) <- c("Pathway", "cluster")
merged <- merge(t, m, all.x=T)

write.csv(merged, file=args[2])


