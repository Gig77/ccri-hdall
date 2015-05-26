options(warn=1)

m <- read.delim("/mnt/projects/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv", stringsAsFactor=F)
m <- m[m$status != "REJECT" & m$non_silent==1 & m$sample != "rem_rel2" & m$var_type == "snp" & m$patient != "E",]
m <- m[m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3"),]
m$sample[m$sample=="rem_dia"] <- "diagnosis"
m$sample[m$sample=="rem_rel"] <- "relapse"

h <- read.delim("/mnt/projects/hdall/results/ras-heterogeneity/ras.hotspots.tsv", stringsAsFactor=F)
h <- h[h$cohort != "non-relapsing",]

merged <- merge(m, h, by=c("patient", "sample", "gene", "chr", "pos", "ref", "alt"), all.x=T, all.y=T)

# Ras pathway mutations per patient
t.dia <- table(h$patient[h$sample=="diagnosis"])
t.rel <- table(h$patient[h$sample=="relapse"])
summary(as.numeric(t.dia))
summary(as.numeric(t.rel))

# number of Ras pathway mutations at diagnosis and relapse
length(unique(h$patient[h$sample=="diagnosis"]))
length(unique(h$patient[h$sample=="relapse"]))

# distribution of Ras pathway mutations at diagnosis and relapse
summary(h$frequency[h$sample=="diagnosis"])
summary(h$frequency[h$sample=="relapse"])

# Ras pathway mutations not detected by MuTect
sum(is.na(merged$impact[merged$sample=="diagnosis"]))
sum(is.na(merged$impact[merged$sample=="relapse"]))

# Ras pathway mutations not detected at hotspots
merged[is.na(merged$frequency) & merged$gene != "FLT3",]
