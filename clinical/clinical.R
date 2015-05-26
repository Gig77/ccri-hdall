options(warn=1)
library(gridExtra)
library(vcd)
library(survival)

rm(list=ls())
source("/mnt/projects/hdall/scripts/clinical/test_pairwise.R")

set.seed(22)

patients_exome <- c("314", "399", "430", "446", "460", "545", "715", "786", "792", "818", "842", "A", "B", "C", "D", "X", "Y", "1021247", "592")
patients_rel_only <- c("1017005", "1021865", "1023545", "AD15", "BL16", "BM18", "CA18", "DM1", "FE1", "FS1", "GD18", "GD1", "HJ15", "HJA15", "HL1", "KA17", "KJ17", "KL16", "LB17", "LM18", "MJ1", "ML10", "MV16", "NS18", "PC16", "RT15", "RT16", "SJM16", "SKR1", "SL1", "SLM1", "ST14", "WA1", "ZE13")
patients_dia_only <- c("1004564", "1010661", "1010781", "1019964", "1020076", "1021087", "1023338", "1023616", "1024589", "1026233", "1026662", "B100", "EF7", "FB14", "G44", "HD7")
patients_non_rel <- c("331", "380", "442", "350", "461", "466", "529", "591", "602", "619", "633", "634", "642", "653", "666", "672", "697", "698", "700", "709", "724", "762", "776", "777", "779", "782", "409", "NRD_1", "73", "NRD_2", "NRD_3", "60", "594", "687", "748", "754", "646", "530", "718", "681", "39", "49", "45", "54", "110", "111", "134", "143", "147", "199", "NRD_4")

c <- read.delim("/mnt/projects/hdall/results/clinical/clinical_data.tsv", na.strings=c("", "NA", "n/a", "n/d", " ", "early (CNS)"))
c <- c[!(c$patient_id %in% c("E", "RN14046", "1025678", "1021186")),]
c <- c[!c$source %in% c("REMOVE", "not included"),]

#----
# CLEAN UP DATA
#----

c$age_dia <- as.numeric(c$age_dia)
c$sex[c$sex!="m" & c$sex!="f"] <- NA
c$source <- as.factor(as.character(c$source))
c$sex <- as.factor(as.character(c$sex))
c$blasts_dia <- as.numeric(as.character(c$blasts_dia))
c$blasts_rel <- as.numeric(as.character(c$blasts_rel))
c$first_rem_months <- as.numeric(as.character(c$first_rem_months))
c$second_rem_months <- as.numeric(as.character(c$second_rem_months))
c$age_rel <- ifelse(c$cohort=="non-relapsing", NA, c$age_dia+c$first_rem_months/12)
c$rel_timepoint <- as.character(c$rel_timepoint)
c$rel_timepoint[c$rel_timepoint=="very early"] <- "early"
c$rel_timepoint <- as.factor(c$rel_timepoint)
#----
# MERGE MUTATION DATA
#----

m <- read.delim("/mnt/projects/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv")
m <- m[m$status!="REJECT" & m$non_silent==T & m$freq_leu >= 0.05,]

m.merged <- merge(m[m$sample=="rem_dia", c("patient", "chr", "pos", "ref", "alt", "freq_leu", "gene")], m[m$sample=="rem_rel",c("patient", "chr", "pos", "ref", "alt", "freq_leu", "gene")], by=c("patient", "chr", "pos", "ref", "alt", "gene"), all=T, suffixes=c(".dia", ".rel")) 

m.exome <- read.delim("/mnt/projects/hdall/results/filtered-variants.cosmic.normaf.tsv")
m.exome <- m.exome[m.exome$status!="REJECT" & m.exome$freq_leu >= 0.1,]
m.exome.ns <- m.exome[m.exome$non_silent==T,]

m.nonrel <- read.delim("/mnt/projects/hdall/results/reseq/filtered-variants.reseq.nonrel.tsv")
m.nonrel <- m.nonrel[m.nonrel$status!="REJECT" & m.nonrel$non_silent==T & m.nonrel$freq_leu >= 5,]
m.nonrel$freq_leu = m.nonrel$freq_leu / 100 

# TODO: determine number of RAS pathway hotspot mutations from own mutation caller plus FLT3 mutations from MuTect
m.ras <- read.delim("/mnt/projects/hdall/results/ras-heterogeneity/ras.hotspots.tsv", stringsAsFactors=F)

m.ras.flt3 <- m[m$sample=="rem_dia" & m$gene=="FLT3", c("patient", "sample", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
m.ras.flt3[sapply(m.ras.flt3, is.factor)] <- lapply(m.ras.flt3[sapply(m.ras.flt3, is.factor)], as.character) # convert all factors in dataframe to characters
m.ras.flt3$sample[m.ras.flt3$sample=="rem_dia"] <- "diagnosis"
m.ras.flt3$sample[m.ras.flt3$sample=="rem_rel" | m.ras.flt3$sample=="rem_rel2"] <- "relapse"
m.ras.flt3 <- data.frame(patient=m.ras.flt3[,1], cohort="relapsing", m.ras.flt3[,2:ncol(m.ras.flt3)], stringsAsFactors=F)
names(m.ras.flt3) <- names(m.ras)
m.ras <- rbind(m.ras, m.ras.flt3)

m.ras.flt3.nr <- m.nonrel[m.nonrel$gene=="FLT3", c("patient", "gene", "chr", "pos", "ref", "alt", "dp_leu_var", "dp_leu_tot", "freq_leu")]
m.ras.flt3.nr[sapply(m.ras.flt3.nr, is.factor)] <- lapply(m.ras.flt3.nr[sapply(m.ras.flt3.nr, is.factor)], as.character) # convert all factors in dataframe to characters
m.ras.flt3.nr <- data.frame(patient=m.ras.flt3.nr[,1], cohort="non-relapsing", sample="diagnosis", m.ras.flt3.nr[,2:ncol(m.ras.flt3.nr)], stringsAsFactors=F)
names(m.ras.flt3) <- names(m.ras)
m.ras <- rbind(m.ras, m.ras.flt3)

m.ras <- unique(m.ras)  # some FLT3 mutations might now be duplicated

num.ras.dia <- aggregate(frequency~patient, data=m.ras[m.ras$sample=="diagnosis",], FUN=length)
names(num.ras.dia) <- c("patient_id", "num.ras.dia")
num.ras.rel <- aggregate(frequency~patient, data=m.ras[m.ras$sample=="relapse",], FUN=length)
names(num.ras.rel) <- c("patient_id", "num.ras.rel")


#----
# ADD ATTRIBUTES TO CLINICAL TABLE
#----

c$sec_rel <- as.logical(ifelse(c$cohort=="non-relapsing", NA, ifelse(c$second_event_after_first_relapse=="death", NA, c$second_event_after_first_relapse=="2nd rel")))  # if someone died, we don't know whether s/he would have gotten a second relapse; therefore NA
c$had_second_event_after_first <- as.logical(ifelse(c$cohort=="non-relapsing", NA, c$second_event_after_first_relapse=="2nd rel" | c$second_event_after_first_relapse=="death"))
c$dead <- c$overall_outcome == "death"


c$crebbp <- as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="CREBBP", "patient"])
c$crebbp.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$crebbp)
c$crebbp.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP" & m$sample=="rem_dia", "patient"])))
c$crebbp.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
c$crebbp.nonrelapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, c$crebbp, NA)
c$crebbp.dia <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, ifelse((!is.na(c$crebbp.nonrelapsing) & c$crebbp.nonrelapsing) | (!is.na(c$crebbp.relapsing.dia) & c$crebbp.relapsing.dia), TRUE, FALSE))

# allelic frequency CREBBP at diagnosis
af <- aggregate(freq_leu~patient, data=m[m$gene=="CREBBP" & m$sample=="rem_dia", c("patient", "freq_leu")], FUN=max)
names(af) <- c("patient_id", "crebbp.relapsing.dia.af")
c <- merge(c, af, all.x=T)

# merge number of RAS pathway mutations
c <- merge(c, num.ras.dia, all.x=T)
c$num.ras.dia[is.na(c$num.ras.dia)] <- ifelse(as.character(c$patient_id[is.na(c$num.ras.dia)]) %in% c(patients_rel_only), NA, 0)
c <- merge(c, num.ras.rel, all.x=T)
c$num.ras.rel[is.na(c$num.ras.rel)] <- ifelse(as.character(c$patient_id[is.na(c$num.ras.rel)]) %in% c(patients_non_rel, patients_dia_only), NA, 0)

# allelic frequency CREBBP at relapse
af <- aggregate(freq_leu~patient, data=m[m$gene=="CREBBP" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), c("patient", "freq_leu")], FUN=max)
names(af) <- c("patient_id", "crebbp.relapsing.rel.af")
c <- merge(c, af, all.x=T)

c$kras <- as.character(c$patient_id) %in% as.character(m[m$gene=="KRAS", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="KRAS", "patient"])
c$kras.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$kras)
c$kras.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="KRAS" & m$sample=="rem_dia", "patient"])))
c$kras.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="KRAS" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
c$kras.nonrelapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, c$kras, NA)
c$kras.dia <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, ifelse((!is.na(c$kras.nonrelapsing) & c$kras.nonrelapsing) | (!is.na(c$kras.relapsing.dia) & c$kras.relapsing.dia), TRUE, FALSE))

c$nras <- as.character(c$patient_id) %in% as.character(m[m$gene=="NRAS", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="NRAS", "patient"])
c$nras.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$nras)
c$nras.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="NRAS" & m$sample=="rem_dia", "patient"])))
c$nras.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="NRAS" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
c$nras.nonrelapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, c$nras, NA)
c$nras.dia <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, ifelse((!is.na(c$nras.nonrelapsing) & c$nras.nonrelapsing) | (!is.na(c$nras.relapsing.dia) & c$nras.relapsing.dia), TRUE, FALSE))

c$ptpn11 <- as.character(c$patient_id) %in% as.character(m[m$gene=="PTPN11", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="PTPN11", "patient"])
c$ptpn11.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$ptpn11)
c$ptpn11.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="PTPN11" & m$sample=="rem_dia", "patient"])))
c$ptpn11.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="PTPN11" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
c$ptpn11.nonrelapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, c$ptpn11, NA)
c$ptpn11.dia <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, ifelse((!is.na(c$ptpn11.nonrelapsing) & c$ptpn11.nonrelapsing) | (!is.na(c$ptpn11.relapsing.dia) & c$ptpn11.relapsing.dia), TRUE, FALSE))

c$kras.or.ptpn11 <- c$kras | c$ptpn11
c$kras.or.ptpn11.relapsing <- c$kras.relapsing | c$ptpn11.relapsing
c$kras.or.ptpn11.relapsing.dia <- c$kras.relapsing.dia | c$ptpn11.relapsing.dia
c$kras.or.ptpn11.relapsing.rel <- c$kras.relapsing.rel | c$ptpn11.relapsing.rel

#c$crebbp.dia <- ifelse(as.character(c$patient_id) %in% patients_rel_only, NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP" & m$sample=="rem_dia", "patient"])) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="CREBBP", "patient"]))
#c$crebbp.rel <- ifelse(as.character(c$patient_id) %in% patients_dia_only, NA, as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]))
#c$crebbp.all <- as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP", "patient"])

c$ras <- as.character(c$patient_id) %in% as.character(m[m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3"), "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3"), "patient"])
c$ras.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$ras)
c$ras.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3") & m$sample=="rem_dia", "patient"])))
c$ras.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3") & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
c$ras.nonrelapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, c$ras, NA)
c$ras.dia <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, ifelse((!is.na(c$ras.nonrelapsing) & c$ras.nonrelapsing) | (!is.na(c$ras.relapsing.dia) & c$ras.relapsing.dia), TRUE, FALSE))

c$ras.heterogeneous.dia <- ifelse(c$num.ras.dia == 0, NA, c$num.ras.dia > 1)
c$ras.heterogeneous.rel <- ifelse(c$num.ras.rel == 0, NA, c$num.ras.rel > 1)
c$ras.conserved <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only, patients_rel_only, as.character(c[!c$ras,"patient_id"])), NA, (as.character(c$patient_id) %in% as.character(m.merged[m.merged$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3") & !is.na(m.merged$freq_leu.dia) & !is.na(m.merged$freq_leu.rel), "patient"])))


# allelic frequency RAS pathway mutation at diagnosis
af <- aggregate(frequency~patient, data=m.ras[m.ras$sample=="diagnosis", c("patient", "frequency")], FUN=max)
names(af) <- c("patient_id", "ras.dia.af")
c <- merge(c, af, all.x=T)

# allelic frequency RAS pathway mutation at relapse
af <- aggregate(frequency~patient, data=m.ras[m.ras$sample=="relapse", c("patient", "frequency")], FUN=max)
names(af) <- c("patient_id", "ras.rel.af")
c <- merge(c, af, all.x=T)

c$crebbp.and.kras.relapsing.dia <- c$crebbp.relapsing.dia & c$kras.relapsing.dia
c$crebbp.and.kras.relapsing.rel <- c$crebbp.relapsing.rel & c$kras.relapsing.rel
c$crebbp.and.ras.relapsing.dia <- c$crebbp.relapsing.dia & c$ras.relapsing.dia
c$crebbp.and.ras.relapsing.rel <- c$crebbp.relapsing.rel & c$ras.relapsing.rel

c$num.mut.dia <- ifelse(as.character(c$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome$patient==as.character(x) & m.exome$sample=="rem_dia") }), NA)
c$num.mut.rel <- ifelse(as.character(c$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome$patient==as.character(x) & m.exome$sample=="rem_rel") }), NA)
c$num.mut.dia.ns <- ifelse(as.character(c$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome.ns$patient==as.character(x) & m.exome.ns$sample=="rem_dia") }), NA)
c$num.mut.rel.ns <- ifelse(as.character(c$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome.ns$patient==as.character(x) & m.exome.ns$sample=="rem_rel") }), NA)
c$num.mut.rel.excl.patA <- ifelse(c$patient_id == "A", NA, c$num.mut.rel)
c$num.mut.rel.ns.excl.patA <- ifelse(c$patient_id == "A", NA, c$num.mut.rel.ns)

# dominant RAS pw mutation per timepoint
m.ras.dia <- m[m$sample=="rem_dia" & m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3"),]
suppressWarnings(m.ras.dia <- m.ras.dia[ave(m.ras.dia$freq_leu, m.ras.dia$patient, FUN=max)==m.ras.dia$freq_leu, c("patient", "gene")])
names(m.ras.dia) <- c("patient_id", "ras.pw.gene.dia")
c <- merge(c, m.ras.dia, all.x=T)
c$ras.pw.gene.dia <- as.character(c$ras.pw.gene.dia)
c$ras.pw.gene.dia[is.na(c$ras.pw.gene.dia) & !as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only)] <- "WT"
c$ras.pw.gene.dia <- as.factor(c$ras.pw.gene.dia)

m.ras.rel <- m[(m$sample=="rem_rel" | m$sample=="rem_rel2") & m$gene %in% c("KRAS", "NRAS", "PTPN11", "FLT3"),]
suppressWarnings(m.ras.rel <- m.ras.rel[ave(m.ras.rel$freq_leu, m.ras.rel$patient, FUN=max)==m.ras.rel$freq_leu, c("patient", "gene")])
names(m.ras.rel) <- c("patient_id", "ras.pw.gene.rel")
c <- merge(c, m.ras.rel, all.x=T)
c$ras.pw.gene.rel <- as.character(c$ras.pw.gene.rel)
c$ras.pw.gene.rel[is.na(c$ras.pw.gene.rel) & !as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only)] <- "WT"
c$ras.pw.gene.rel <- as.factor(c$ras.pw.gene.rel)

# ALTERNATIVELY, USE DATA FROM HOTSPOT CALLER; BUT: CHANGES KAPLAN-MAIER PLOT FOR PTPN11 positive case
#m.ras.dia <- m.ras[m.ras$sample=="diagnosis",]
#suppressWarnings(m.ras.dia <- m.ras.dia[ave(m.ras.dia$frequency, m.ras.dia$patient, FUN=max)==m.ras.dia$frequency, c("patient", "gene")])
#names(m.ras.dia) <- c("patient_id", "ras.pw.gene.dia")
#c <- merge(c, m.ras.dia, all.x=T)
#c$ras.pw.gene.dia <- as.character(c$ras.pw.gene.dia)
#c$ras.pw.gene.dia[is.na(c$ras.pw.gene.dia) & !as.character(c$patient_id) %in% c(patients_rel_only)] <- "WT"
#c$ras.pw.gene.dia <- as.factor(c$ras.pw.gene.dia)

#m.ras.rel <- m.ras[m.ras$sample=="relapse",]
#suppressWarnings(m.ras.rel <- m.ras.rel[ave(m.ras.rel$frequency, m.ras.rel$patient, FUN=max)==m.ras.rel$frequency, c("patient", "gene")])
#names(m.ras.rel) <- c("patient_id", "ras.pw.gene.rel")
#c <- merge(c, m.ras.rel, all.x=T)
#c$ras.pw.gene.rel <- as.character(c$ras.pw.gene.rel)
#c$ras.pw.gene.rel[is.na(c$ras.pw.gene.rel) & !as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only)] <- "WT"
#c$ras.pw.gene.rel <- as.factor(c$ras.pw.gene.rel)

#----
# GENE DELETIONS
#----

# data noisy or blast cound too high
inconclusive.patients.dia <- c("1020076", "1023616", "KA14651")
inconclusive.patients.rel <- c("1021631", "1022914", "1025409", "BM18", "FS1", "LM18", "PC16", "SJM16")

c$ikzf1.del.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only, inconclusive.patients.dia), NA, ifelse(c$patient_id %in% c("KE12025", "HD7", "B100", "243"), TRUE, FALSE))
c$ikzf1.del.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only, inconclusive.patients.rel), NA, ifelse(c$patient_id %in% c("446", "MB2", "KE12025", "HJ15", "KA17", "KL16"), TRUE, FALSE))
c$ikzf1.del <- c$ikzf1.del.dia | c$ikzf1.del.rel
#c$ikzf2.del.rel <- ifelse(c$cohort=="relapsing" & !(c$patient_id %in% inconclusive.patients.rel), ifelse(c$patient_id %in% c("Y", "LB17", "KD20493", "CA18", "G", "X", "BL16"), TRUE, FALSE), NA)
#c$ikzf.del.rel <- c$ikzf1.del.rel | c$ikzf2.del.rel
#c$ikzf.del <- c$ikzf1.del | c$ikzf2.del.rel
c$crebbp.relapsing.rel.mut_or_del <- ifelse(is.na(c$crebbp.relapsing.rel), NA, c$crebbp.relapsing.rel | c$patient_id %in% c("C", "LB17", "KE12025", "1024589"))
c$crebbp.relapsing.dia.mut_or_del <- ifelse(is.na(c$crebbp.relapsing.dia), NA, c$crebbp.relapsing.dia | c$patient_id %in% c("592", "KE12025", "HD7"))
c$crebbp.relapsing.mut_or_del <- c$crebbp.relapsing.dia.mut_or_del | c$crebbp.relapsing.rel.mut_or_del

c$kras.or.ikzf1.rel <- c$kras.relapsing.rel | c$ikzf1.del.rel
c$kras.or.nras.or.ikzf1.rel <- c$kras.relapsing.rel | c$nras.relapsing.rel | c$ikzf1.del.rel

cr <- c[c$cohort=="relapsing",]

#----
# CLONAL KINETICS
#----
c$clonal.kinetic <- NA
c$clonal.kinetic[c$patient_id %in% c("430", "446", "818", "B", "C", "842", "D", "1021247", "792", "A")] <- "progression"
c$clonal.kinetic[c$patient_id %in% c("460", "545", "715", "592", "786", "X", "Y")] <- "ancestral-derived"
c$clonal.kinetic <- as.factor(c$clonal.kinetic)

#boxplot(first_rem_months ~ sex, data=cr, na.action=na.exclude, outline=F, names=c("female", "male"))
#stripchart(first_rem_months ~ sex, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)

write.table(c, file="/mnt/projects/hdall/results/clinical/clinical_data.processed.tsv", col.names=T, row.names=F, sep="\t", quote=F)

stop("OK")

#---
# associations with cohort (relapsing vs. non-relapsing), all data
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-cohort-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, include=c("cohort"), exclude=c("patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"))
dev.off()

#---
# associations of cohort (relapsing vs. non-relapsing), mutation data
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-cohort-vs-mutation.pdf")
test_pairwise_assoc(c, include.group=c("cohort", "crebbp.dia", "kras.dia", "nras.dia", "ptpn11.dia", "ras.dia"), exclude.group=list(c("crebbp.dia", "kras.dia", "nras.dia", "ptpn11.dia", "ras.dia")))
dev.off()

#---
# associations of mutation kinetics (progression vs. ancestral clone)
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-kinetic-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, include=c("clonal.kinetic"), exclude=c("patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"))
dev.off()

#---
# KM plots with mutation status at diagnosis
#---


# AT POINT OF DIAGNOSIS, IT IS CURRENTLY NOT CLEAR WHAT WE WANT TO SHOW WITH KM PLOTS (TOTAL REMISSION TIME, FIRST REMISSION TIME???)

#c$relapsing <- c$cohort=="relapsing"

#pdf("/mnt/projects/hdall/results/clinical/kaplan-mutations-at-diagnosis.pdf")

#plot(survfit(Surv(time=first_rem_months, relapsing)~kras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pOS", conf.int=F)
#legend(170, 0.2, c("wtKRAS dia", "mKRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, had_second_event_after_first)~kras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pEFS", conf.int=F)
#legend(170, 0.2, c("wtKRAS dia", "mKRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, relapsing)~nras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pOS", conf.int=F)
#legend(170, 0.2, c("wtNRAS dia", "mNRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, had_second_event_after_first)~nras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pEFS", conf.int=F)
#legend(170, 0.2, c("wtNRAS dia", "mNRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, relapsing)~ptpn11.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pOS", conf.int=F)
#legend(170, 0.2, c("wtPTPN11 dia", "mPTPN11 dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, had_second_event_after_first)~ptpn11.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pEFS", conf.int=F)
#legend(170, 0.2, c("wtPTPN11 dia", "mPTPN11 dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, relapsing)~ras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pOS", conf.int=F)
#legend(170, 0.2, c("wtRAS dia", "mRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, had_second_event_after_first)~ras.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pEFS", conf.int=F)
#legend(170, 0.2, c("wtRAS dia", "mRAS dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, relapsing)~crebbp.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pOS", conf.int=F)
#legend(170, 0.2, c("wtCREBBP dia", "mCREBBP dia"), lwd=c(1,1), col=c("blue", "red"))

#plot(survfit(Surv(time=first_rem_months, had_second_event_after_first)~crebbp.dia, data=c), col=c("blue", "red"), xlab="remission (months)", ylab="pEFS", conf.int=F)
#legend(170, 0.2, c("wtCREBBP dia", "mCREBBP dia"), lwd=c(1,1), col=c("blue", "red"))

#dev.off()

#---
# associations with CREBBP mutation status at diagnosis
#---
c$crebbp.dia.deleterious <- ifelse(as.character(c$patient_id) %in% c(patients_rel_only), NA, as.character(c$patient_id) %in% as.character(m[m$gene=="CREBBP" & m$sample=="rem_dia" & m$deleterious=="yes", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="CREBBP" & m.nonrel$deleterious=="yes", "patient"]))

pdf("/mnt/projects/hdall/results/clinical/clinical-crebbp-dia-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, 
		include=c("crebbp.dia", "crebbp.relapsing.dia", "crebbp.dia.deleterious"), 
		exclude=c("crebbp", "crebbp.relapsing", "crebbp.relapsing.rel", "crebbp.nonrelapsing", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia", "crebbp.and.kras.relapsing.rel", "crebbp.relapsing.mut_or_del", "crebbp.relapsing.dia.mut_or_del", "crebbp.relapsing.rel.mut_or_del", 
				  "patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"),
		exclude.group=list(c("crebbp.dia", "crebbp.relapsing.dia", "crebbp.dia.deleterious")))
dev.off()

#---
# associations with CREBBP mutation status at relapse
#---
c$crebbp.relapsing.rel.no_dis <- c$crebbp.relapsing.rel 
c$crebbp.relapsing.rel.no_dis[c$patient_id %in% patients_exome] <- NA

pdf("/mnt/projects/hdall/results/clinical/clinical-crebbp-rel-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, 
		include=c("crebbp.relapsing.rel", "crebbp.relapsing.rel.no_dis"), 
		exclude=c("crebbp", "crebbp.relapsing", "crebbp.dia", "crebbp.relapsing.dia", "crebbp.nonrelapsing", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia", "crebbp.and.kras.relapsing.rel", "crebbp.relapsing.mut_or_del", "crebbp.relapsing.dia.mut_or_del", "crebbp.relapsing.rel.mut_or_del", 
				  "patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"),
	    exclude.group=list(c("wtkras.wtcrebbp.relapsing.rel", "crebbp.relapsing.rel", "wtkras.mcrebbp.relapsing.rel", "mkras.mcrebbp.relapsing.rel", "mkras.wtcrebbp.relapsing.rel", "crebbp.relapsing.rel.no_dis"))  
)
dev.off()

#---
# associations with CREBBP+/KRAS+ mutation status at relapse
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-crebbp+kras-rel-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, 
		include=c("crebbp.and.kras.relapsing.rel"), 
		exclude=c("crebbp", "crebbp.relapsing", "crebbp.dia", "crebbp.relapsing.dia", "crebbp.relapsing.rel", "crebbp.nonrelapsing", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia", "crebbp.relapsing.mut_or_del", "crebbp.relapsing.dia.mut_or_del", "crebbp.relapsing.rel.mut_or_del", "kras", 
				  "kras.relapsing", "kras.relapsing.dia", "kras.relapsing.rel", "kras.nonrelapsing", "kras.dia", "ras", "ras.relapsing", "ras.relapsing.dia", "ras.relapsing.rel", "ras.nonrelapsing", "ras.dia", "num.ras.rel", "num.ras.dia", "kras.or.ptpn11", "kras.or.ptpn11.relapsing", "kras.or.ptpn11.relapsing.dia", "kras.or.ptpn11.relapsing.rel", "kras.or.ikzf1.rel", "kras.or.nras.or.ikzf1.rel", 
				  "patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment")
)
dev.off()

#---
# associations with RAS pw genes at diagnosis and relapse
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-ras.pw.gene-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, 
		include=c("ras.pw.gene.dia", "ras.pw.gene.rel"), 
		exclude=c("patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"),
		exclude.group=list(c("ras.pw.gene.dia", "ras.pw.gene.rel", "kras.relapsing", "kras.relapsing.dia", "kras.relapsing.rel", "kras.nonrelapsing", "kras.dia", "ras", "ras.relapsing", "ras.relapsing.dia", "ras.relapsing.rel", "ras.nonrelapsing", "ras.dia", "num.ras.rel", "num.ras.dia", "kras.or.ptpn11", "kras.or.ptpn11.relapsing", "kras.or.ptpn11.relapsing.dia", "kras.or.ptpn11.relapsing.rel", "kras.or.ikzf1.rel", "kras.or.nras.or.ikzf1.rel", "nras", "nras.relapsing", "nras.relapsing.rel", "nras.relapsing.dia", "nras.dia", "ptpn11", "ptpn11.dia", "ptpn11.relapsing", "ptpn11.relapsing.rel", "ptpn11.relapsing.dia", "ptpn11.dia", "kras", "nras", "nras.relapsing", "ptpn11", "ptpn11.relapsing", "crebbp.and.kras.relapsing.rel", "kras", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia")))
dev.off()

#---
# associations with PTPN11 mutation status
#---
pdf("/mnt/projects/hdall/results/clinical/clinical-ptpn11-vs-all.pdf")
test_pairwise_assoc(c, sig.level=0.99, 
		include=c("ptpn11.relapsing.rel", "ptpn11.relapsing.dia", "ptpn11.dia"), 
		exclude=c("ptpn11", "ptpn11.relapsing", "ptpn11.nonrelapsing", "ras.pw.gene.rel", "ras.pw.gene.dia", "kras.or.ptpn11.relapsing.rel", "kras.or.ptpn11", "kras.or.ptpn11.relapsing", "ras.relapsing.rel", "ras.relapsing.dia", "kras.or.ptpn11.relapsing.dia", "ras", "ras.relapsing", "ras.dia",
				"patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"),  
		exclude.group=list(c("ptpn11.relapsing.rel", "ptpn11.relapsing.dia", "ptpn11.dia"))
)
dev.off()

#---
# association of RAS pw clonal status at diagnosis with cohort
#---
pdf("/mnt/projects/hdall/results/clinical/ras.pw.af.dia-relapsing.vs.nonrelapsing.pdf")
par(mfrow=c(2,3))
boxplot(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="KRAS",], xlab=sprintf("p=%.2g", kruskal.test(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="KRAS",], na.action=na.exclude)$p.value), ylab="AF", main="KRAS", na.action=na.exclude, outline=F, cex.axis=0.8)
stripchart(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="KRAS",], method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c$ras.pw.gene.dia, add=T)
boxplot(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="NRAS",], xlab=sprintf("p=%.2g", kruskal.test(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="NRAS",], na.action=na.exclude)$p.value), ylab="AF", main="NRAS", na.action=na.exclude, outline=F, cex.axis=0.8)
stripchart(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="NRAS",], method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c$ras.pw.gene.dia, add=T)
boxplot(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="PTPN11",], xlab=sprintf("p=%.2g", kruskal.test(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="PTPN11",], na.action=na.exclude)$p.value), ylab="AF", main="PTPN11", na.action=na.exclude, outline=F, cex.axis=0.8)
stripchart(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="PTPN11",], method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c$ras.pw.gene.dia, add=T)
boxplot(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="FLT3",], xlab=sprintf("p=%.2g", kruskal.test(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="FLT3",], na.action=na.exclude)$p.value), ylab="AF", main="FLT3", na.action=na.exclude, outline=F, cex.axis=0.8)
stripchart(ras.dia.af~cohort, c[c$ras.pw.gene.dia=="FLT3",], method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c$ras.pw.gene.dia, add=T)
boxplot(ras.dia.af~cohort, c, xlab=sprintf("p=%.2g", kruskal.test(ras.dia.af~cohort, c, na.action=na.exclude)$p.value), ylab="AF", main="combined", na.action=na.exclude, outline=F, cex.axis=0.8)
stripchart(ras.dia.af~cohort, c, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c$ras.pw.gene.dia, add=T)
dev.off()

# correlation of mutations with age
#---
png("/mnt/projects/hdall/results/clinical/supp-figure.corr.age-dia.num-mut-dia.png", width=2048, height=1024)
par(mfrow=c(1,2), cex=2)

fit <- lm(num.mut.dia~age_dia, data=c)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(c$age_dia, c$num.mut.dia, xlab="age at diagnosis (years)", ylab="Number of mutations", xlim=c(0, 20), ylim=c(0, 250), main=sprintf("R=%.2f, p=%.2g", R, p))
abline(fit, col="red")

fit <- lm(num.mut.rel.excl.patA~age_rel, data=c)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(c$age_rel, c$num.mut.rel.excl.patA, xlab="age at relapse (years)", ylab="Number of mutations", xlim=c(0, 20), ylim=c(0, 250), main=sprintf("R=%.2f, p=%.2g", R, p))
abline(fit, col="red")

dev.off()

pdf("/mnt/projects/hdall/results/clinical/clinical.pdf")
tests <- test_pairwise_assoc(c, 
		sig.level=0.1, 
		exclude=c("patient_id", "exome", "panel", "mrd_level_rel", "rel_protocol", "BM.transplantation.date", "study_no_relapse", "patno_kiel_rem", "patno_berlin_rel", "X", "comment"),
		exclude.group=list(c("crebbp", "crebbp.relapsing", "crebbp.relapsing.dia", "crebbp.relapsing.rel", "crebbp.nonrelapsing", "crebbp.dia", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia", "crebbp.and.kras.relapsing.rel", "crebbp.relapsing.mut_or_del", "crebbp.relapsing.dia.mut_or_del", "crebbp.relapsing.rel.mut_or_del"),
				 		   c("kras", "kras.relapsing", "kras.relapsing.dia", "kras.relapsing.rel", "kras.nonrelapsing", "kras.dia", "ras", "ras.relapsing", "ras.relapsing.dia", "ras.relapsing.rel", "ras.nonrelapsing", "ras.dia", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "crebbp.and.kras.relapsing.dia", "crebbp.and.kras.relapsing.rel", "num.ras.rel", "num.ras.dia", "kras.or.ptpn11", "kras.or.ptpn11.relapsing", "kras.or.ptpn11.relapsing.dia", "kras.or.ptpn11.relapsing.rel", "kras.or.ikzf1.rel", "kras.or.nras.or.ikzf1.rel"),
						   c("nras", "nras.relapsing", "nras.relapsing.dia", "nras.relapsing.rel", "nras.nonrelapsing", "nras.dia", "ras", "ras.relapsing", "ras.relapsing.dia", "ras.relapsing.rel", "ras.nonrelapsing", "ras.dia", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "num.ras.rel", "num.ras.dia", "kras.or.nras.or.ikzf1.rel"),
						   c("ptpn11", "ptpn11.relapsing", "ptpn11.relapsing.dia", "ptpn11.relapsing.rel", "ptpn11.nonrelapsing", "ptpn11.dia", "ras", "ras.relapsing", "ras.relapsing.dia", "ras.relapsing.rel", "ras.nonrelapsing", "ras.dia", "crebbp.and.ras.relapsing.dia", "crebbp.and.ras.relapsing.rel", "num.ras.rel", "num.ras.dia", "kras.or.ptpn11", "kras.or.ptpn11.relapsing", "kras.or.ptpn11.relapsing.dia", "kras.or.ptpn11.relapsing.rel"),
						   c("num.mut.dia", "num.mut.dia.ns"),
						   c("ikzf1.del.dia", "ikzf1.del.rel", "ikzf1.del", "kras.or.ikzf1.rel", "kras.or.nras.or.ikzf1.rel"),
						   c("kras.or.ikzf1.rel", "kras.or.nras.or.ikzf1.rel"),
						   c("num.mut.rel", "num.mut.rel.ns", "num.mut.rel.excl.patA", "num.mut.rel.ns.excl.patA"),
						   c("second_event_after_first_relapse", "overall_outcome", "dead", "sec_rel", "had_second_event_after_first"))
)
dev.off()


#----------------------------------------------------------

c$mkras.wtcrebbp.relapsing.rel <- c$kras.relapsing.rel & !c$crebbp.relapsing.rel
c$wtkras.mcrebbp.relapsing.rel <- !c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$mkras.mcrebbp.relapsing.rel <- c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$wtkras.wtcrebbp.relapsing.rel <- !c$kras.relapsing.rel & !c$crebbp.relapsing.rel

pdf("/mnt/projects/hdall/results/clinical/ras-crebbp-mrd.pdf")

pair.table <- table(c[,c("mrd_risk_rel", "mkras.wtcrebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "wtkras.mcrebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "mkras.mcrebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "wtkras.wtcrebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "ras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "kras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "nras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "ptpn11.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

c$flt3 <- as.character(c$patient_id) %in% as.character(m[m$gene=="FLT3", "patient"]) | as.character(c$patient_id) %in% as.character(m.nonrel[m.nonrel$gene=="FLT3", "patient"])
c$flt3.relapsing <- ifelse(as.character(c$patient_id) %in% patients_non_rel, NA, c$flt3)
c$flt3.relapsing.dia <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_rel_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="FLT3" & m$sample=="rem_dia", "patient"])))
c$flt3.relapsing.rel <- ifelse(as.character(c$patient_id) %in% c(patients_non_rel, patients_dia_only), NA, (as.character(c$patient_id) %in% as.character(m[m$gene=="FLT3" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"])))
pair.table <- table(c[,c("mrd_risk_rel", "flt3.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("mrd_risk_rel", "crebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

dev.off()

#----------------------------------------------------------
pdf("/mnt/projects/hdall/results/clinical/ras-crebbp.time-to-relapse.pdf")

pair.table <- table(c[,c("rel_timepoint", "ras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "kras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "nras.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "ptpn11.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "flt3.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "crebbp.relapsing.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

dev.off()



# association of CREBBP/KRAS double mutants at relapse in a mCREBBP and mKRAS background, respectively
#---
pdf("/mnt/projects/hdall/results/clinical/ras-crebbp-double.vs.single.pdf")

c$mkras_in_mcrebb_background.rel <- ifelse(c$crebbp.relapsing.rel, ifelse(c$kras.relapsing.rel, TRUE, FALSE), NA)

pair.table <- table(c[,c("mrd_risk_rel", "mkras_in_mcrebb_background.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "mkras_in_mcrebb_background.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

c$mcrebbp_in_mkras_background.rel <- ifelse(c$kras.relapsing.rel, ifelse(c$crebbp.relapsing.rel, TRUE, FALSE), NA)

pair.table <- table(c[,c("mrd_risk_rel", "mcrebbp_in_mkras_background.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)

pair.table <- table(c[,c("rel_timepoint", "mcrebbp_in_mkras_background.rel")])
mosaic(pair.table, pop=F, main=sprintf("p=%.2g", fisher.test(pair.table)$p.value))
labeling_cells(text=pair.table)(pair.table)


dev.off()

stop("DONE")

# KAPLAN MAYER CURVES
# plot(survfit(Surv(first_rem_months)~sex, data=c), col=c("blue", "red"), xlab="1st remission (months)", ylab="sex", conf.int=F)
# legend(80, 1, c("female", "male"), lwd=c(1,1), col=c("blue", "red"))
# p <- survdiff(Surv(first_rem_months)~sex, data=c))  # don't know how this works yet
#par(mfrow=c(3,3), mar=c(2.5,4,1.5,2))

#----
# RAS PW AF AT DIA WITH RELAPSE TP
#----

pdf("/mnt/projects/hdall/results/clinical/ras-af-dia.vs.rel-timepoint.pdf")
test_pairwise_assoc(c, include.group=c("ras.dia.af", "rel_timepoint"))
dev.off()

#----
# RAS PW AF AT DIA WITH RELAPSE TP
#----

pdf("/mnt/projects/hdall/results/clinical/ras-heterogeneity.pdf")
test_pairwise_assoc(c[c$cohort=="relapsing",], include.group=c("ras.heterogeneous.dia", "rel_timepoint"))
test_pairwise_assoc(c[c$cohort=="relapsing",], include.group=c("ras.heterogeneous.dia", "first_rem_months"))
test_pairwise_assoc(c[c$cohort=="relapsing",], include.group=c("ras.heterogeneous.dia", "mrd_risk_dia"))
test_pairwise_assoc(c[c$cohort=="relapsing",], include.group=c("ras.heterogeneous.dia", "mrd_risk_rel"))
test_pairwise_assoc(c, include.group=c("ras.heterogeneous.rel", "rel_timepoint"))
test_pairwise_assoc(c, include.group=c("ras.heterogeneous.rel", "first_rem_months"))
test_pairwise_assoc(c, include.group=c("ras.heterogeneous.rel", "second_rem_months"))
test_pairwise_assoc(c, include.group=c("ras.heterogeneous.rel", "mrd_risk_dia"))
test_pairwise_assoc(c, include.group=c("ras.heterogeneous.rel", "mrd_risk_rel"))
dev.off()

#----
# FIRST REMISSION DURATION BY CREBBP MUTATION STATUS
#----

boxplot(first_rem_months ~ crebbp.all, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP", "mCREBBP"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ crebbp.all, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$crebbp.all==T], cr$first_rem_months[cr$crebbp.all==F])$p.value)))

boxplot(first_rem_months ~ crebbp.dia, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP dia", "mCREBBP dia"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ crebbp.dia, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$crebbp.dia==T], cr$first_rem_months[cr$crebbp.dia==F])$p.value)))

boxplot(first_rem_months ~ crebbp.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP rel", "mCREBBP rel"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ crebbp.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$crebbp.rel==T], cr$first_rem_months[cr$crebbp.rel==F])$p.value)))

#----
# FIRST REMISSION DURATION BY RAS PW MUTATION
#----

boxplot(first_rem_months ~ ras.all, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS", "mRAS"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ ras.all, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$ras.all==T], cr$first_rem_months[cr$ras.all==F])$p.value)))

boxplot(first_rem_months ~ ras.dia, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS dia", "mRAS dia"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ ras.dia, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$ras.dia==T], cr$first_rem_months[cr$ras.dia==F])$p.value)))

boxplot(first_rem_months ~ ras.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS rel", "mRAS rel"), ylab="1st remission duration (months)", cex.axis=0.8)
stripchart(first_rem_months ~ ras.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$ras.rel==T], cr$first_rem_months[cr$ras.rel==F])$p.value)))

#----
# SECOND REMISSION DURATION BY CREBBP AND RAS MUTATION STATUS
#----

boxplot(second_rem_months ~ crebbp.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP rel", "mCREBBP rel"), ylab="2nd remission duration (months)", cex.axis=0.8)
stripchart(second_rem_months ~ crebbp.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$crebbp.rel==T], cr$first_rem_months[cr$crebbp.rel==F])$p.value)))

boxplot(second_rem_months ~ ras.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS rel", "mRAS rel"), ylab="2nd remission duration (months)", cex.axis=0.8)
stripchart(second_rem_months ~ ras.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 80, paste0("p=", sprintf("%.3f", wilcox.test(cr$first_rem_months[cr$ras.rel==T], cr$first_rem_months[cr$ras.rel==F])$p.value)))

#----
# AGE AT DIAGNOSIS BY CREBBP MUTATION STATUS
#----

boxplot(age_dia ~ crebbp.all, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP", "mCREBBP"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ crebbp.all, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 12, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$crebbp.all==T], cr$age_dia[cr$crebbp.all==F])$p.value)))

boxplot(age_dia ~ crebbp.dia, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP dia", "mCREBBP dia"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ crebbp.dia, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 13, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$crebbp.dia==T], cr$age_dia[cr$crebbp.dia==F])$p.value)))

boxplot(age_dia ~ crebbp.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtCREBBP rel", "mCREBBP rel"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ crebbp.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 9, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$crebbp.rel==T], cr$age_dia[cr$crebbp.rel==F])$p.value)))

#par(mfrow=c(2,3), mar=c(2.5,4,1.5,2))

#----
# AGE AT DIAGNOSIS BY RAS PW MUTATION STATUS
#----

boxplot(age_dia ~ ras.all, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS", "mRAS"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ ras.all, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 12, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$ras.all==T], cr$age_dia[cr$ras.all==F])$p.value)))

boxplot(age_dia ~ ras.dia, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS dia", "mRAS dia"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ ras.dia, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 13, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$ras.dia==T], cr$age_dia[cr$ras.dia==F])$p.value)))

boxplot(age_dia ~ ras.rel, data=cr, na.action=na.exclude, outline=F, names=c("wtRAS rel", "mRAS rel"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ ras.rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 9, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$ras.rel==T], cr$age_dia[cr$ras.rel==F])$p.value)))

#----
# AGE AT DIAGNOSIS BY OUTCOME
#----

boxplot(age_dia ~ death, data=cr, na.action=na.exclude, outline=F, names=c("alive", "dead"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ death, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 15, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$death==T], cr$age_dia[cr$death==F])$p.value)))

boxplot(age_dia ~ sec_rel, data=cr, na.action=na.exclude, outline=F, names=c("no 2nd rel", "with 2nd rel"), ylab="age at diagnosis", cex.axis=0.8)
stripchart(age_dia ~ sec_rel, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)
text(1.5, 15, paste0("p=", sprintf("%.3f", wilcox.test(cr$age_dia[cr$sec_rel==T], cr$age_dia[cr$sec_rel==F])$p.value)))

#----
# OVERALL NUMBER OF MUTATIONS
#----

#par(mfrow=c(2, 2), mar=c(4,4,3,2))

mod <- lm(num.mut.dia~age_dia, data=cr)
plot(cr$age_dia, cr$num.mut.dia, main=sprintf("R=%.2f, p=%.3f", summary(mod)$r.squared, anova(mod)$'Pr(>F)'[1]))
abline(mod, col="red")

mod <- lm(num.mut.rel[patient_id != "A"]~age_rel[patient_id != "A"], data=cr)
plot(cr$age_rel[cr$patient_id != "A"], cr$num.mut.rel[cr$patient_id != "A"], main=sprintf("R=%.2f, p=%.3f", summary(mod)$r.squared, anova(mod)$'Pr(>F)'[1]))
abline(mod, col="red")

mod <- lm(num.mut.rel[patient_id != "A"]~first_rem_months[patient_id != "A"], data=cr)
plot(cr$first_rem_months[cr$patient_id != "A"], cr$num.mut.rel[cr$patient_id != "A"], main=sprintf("R=%.2f, p=%.3f", summary(mod)$r.squared, anova(mod)$'Pr(>F)'[1]))
abline(mod, col="red")

#---
# CATEGORIAL PAIRWISE ASSOCIATIONS (MOSAIC PLOTS)
#---

# par mfrow does not work with mosaic plots!!! argh... one plot per page

mrd.crebbp.relapse <- with(cr, table(mrd_risk_rel, crebbp.rel))
mosaic(mrd.crebbp.relapse, pop=F, sub=sprintf("p=%.3f", fisher.test(mrd.crebbp.relapse)$p.value))
labeling_cells(text=mrd.crebbp.relapse)(mrd.crebbp.relapse)

mrd.ras.relapse <- with(cr, table(mrd_risk_rel, ras.rel))
mosaic(mrd.ras.relapse, pop=F, sub=sprintf("p=%.3f", fisher.test(mrd.ras.relapse)$p.value))
labeling_cells(text=mrd.ras.relapse)(mrd.ras.relapse)

mrd.ras.crebbp.relapse <- with(cr, table(mrd_risk_rel, ras.crebbp.rel))
mosaic(mrd.ras.crebbp.relapse, pop=F, sub=sprintf("p=%.3f", fisher.test(mrd.ras.crebbp.relapse)$p.value))
labeling_cells(text=mrd.ras.crebbp.relapse)(mrd.ras.crebbp.relapse)

death.crebbp <- with(cr, table(death, crebbp.all))
mosaic(death.crebbp, pop=F, sub=sprintf("p=%.3f", fisher.test(death.crebbp)$p.value))
labeling_cells(text=death.crebbp)(death.crebbp)

death.crebbp.rel <- with(cr, table(death, crebbp.rel))
mosaic(death.crebbp.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(death.crebbp.rel)$p.value))
labeling_cells(text=death.crebbp.rel)(death.crebbp.rel)

death.ras <- with(cr, table(death, ras.all))
mosaic(death.ras, pop=F, sub=sprintf("p=%.3f", fisher.test(death.ras)$p.value))
labeling_cells(text=death.ras)(death.ras)

death.ras.rel <- with(cr, table(death, ras.rel))
mosaic(death.ras.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(death.ras.rel)$p.value))
labeling_cells(text=death.ras.rel)(death.ras.rel)

death.ras.crebbp.rel <- with(cr, table(death, ras.crebbp.rel))
mosaic(death.ras.crebbp.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(death.ras.crebbp.rel)$p.value))
labeling_cells(text=death.ras.crebbp.rel)(death.ras.crebbp.rel)

crebbp.kras.dia <- with(cr, table(crebbp.dia, kras.dia))
mosaic(crebbp.kras.dia, pop=F, sub=sprintf("p=%.3f", fisher.test(crebbp.kras.dia)$p.value))
labeling_cells(text=crebbp.kras.dia)(crebbp.kras.dia)

crebbp.kras.rel <- with(cr, table(crebbp.rel, kras.rel))
mosaic(crebbp.kras.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(crebbp.kras.rel)$p.value))
labeling_cells(text=crebbp.kras.rel)(crebbp.kras.rel)

secrel.crebbp.rel <- with(cr, table(crebbp.rel, sec_rel))
mosaic(secrel.crebbp.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(secrel.crebbp.rel)$p.value))
labeling_cells(text=secrel.crebbp.rel)(secrel.crebbp.rel)

secrel.ras.rel <- with(cr, table(ras.rel, sec_rel))
mosaic(secrel.ras.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(secrel.ras.rel)$p.value))
labeling_cells(text=secrel.ras.rel)(secrel.ras.rel)

secrel.ras.crebbp.rel <- with(cr, table(ras.crebbp.rel, sec_rel))
mosaic(secrel.ras.crebbp.rel, pop=F, sub=sprintf("p=%.3f", fisher.test(secrel.ras.crebbp.rel)$p.value))
labeling_cells(text=secrel.ras.crebbp.rel)(secrel.ras.crebbp.rel)

ras.risk.dia <- with(cr, table(kras.dia, mrd_risk_dia))
mosaic(ras.risk.dia, pop=F)
labeling_cells(text=ras.risk.dia)(ras.risk.dia)

dev.off()
