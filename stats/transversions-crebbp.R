min.af <- 0.20
min.dp.leu <- 20
min.dp.rem <- 20

t <- read.delim("/mnt/projects/hdall/results/filtered-variants.tsv", stringsAsFactors=F)
t <- t[t$status != "REJECT" & t$var_type == "snp" & t$freq_leu>=min.af & t$dp_rem_tot>=min.dp.rem & t$dp_leu_tot>=min.dp.leu,]
names(t)[names(t)=="sample"] <- "timepoint"

mcrebbp.dia <- c("545", "C", "Y")
mcrebbp.rel <- c("545", "715", "842", "A", "C", "Y", "399")

t$gt_transversion <- FALSE
t$gt_transversion[(t$ref=="G" & t$alt=="T")|(t$ref=="C" & t$alt=="A")] <- TRUE

mcrebbp.dia <- c("545", "C", "Y")
mcrebbp.rel <- c("545", "715", "842", "A", "C", "Y", "399")

t$mcrebbp <- FALSE
t$mcrebbp[t$timepoint=="rem_dia" & t$patient %in% mcrebbp.dia] <- TRUE
t$mcrebbp[t$timepoint=="rem_rel" & t$patient %in% mcrebbp.rel] <- TRUE

# 1. Is there a statistical difference between the frequency of G-T transversion somatic 
# mutations in diagnostic samples compared to relapse samples in patients who gain a 
# CREBBP mutation at relapse?
		
t.gained <- t[!t$patient %in% mcrebbp.dia & t$patient %in% mcrebbp.rel,]
print(unique(t.gained$patient))
cont <- table(t.gained[,c("timepoint", "gt_transversion")]) ; print(cont)
fisher.test(cont)$p.value

# 2. Is there a statistical difference between the frequency of G-T transversion 
# somatic mutations in diagnostic samples compared to  relapse samples for patients 
# who have CREBBP mutation at diagnosis and relapse?

t.cons <- t[t$patient %in% mcrebbp.dia & t$patient %in% mcrebbp.rel,]
print(unique(t.cons$patient))
cont <- table(t.cons[,c("timepoint", "gt_transversion")]) ; print(cont)
fisher.test(cont)$p.value

# 3. Is there a statistical difference between the frequency of G-T transversion somatic 
# mutations in diagnostic samples with CREBBP mutations compared to diagnostic samples 
# that are CREBBP wildtype?

cont <- table(t[t$timepoint=="rem_dia", c("mcrebbp", "gt_transversion")]) ; print(cont)
fisher.test(cont)$p.value

# 4. Is there a statistical difference between the frequency of G-T transversion somatic 
# mutations  in relapse samples with CREBBP mutations compared to relapse samples that are 
# CREBBP  wildtype?

cont <- table(t[t$timepoint=="rem_rel", c("mcrebbp", "gt_transversion")]) ; print(cont)
fisher.test(cont)$p.value

