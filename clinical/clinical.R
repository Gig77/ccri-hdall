library(gridExtra)

rm(list=ls())
set.seed(22)

patients_exome <- c("314", "399", "430", "446", "460", "545", "715", "786", "792", "818", "842", "A", "B", "C", "D", "X", "Y", "1021247", "592")
patients_rel_only <- c("1017005", "1021865", "1023545", "AD15", "BL16", "BM18", "CA18", "DM1", "FE1", "FS1", "GD18", "GD1", "HJ15", "HJA15", "HL1", "KA17", "KJ17", "KL16", "LB17", "LM18", "MJ1", "ML10", "MV16", "NS18", "PC16", "RT15", "RT16", "SJM16", "SKR1", "SL1", "SLM1", "ST14", "WA1", "ZE13")
patients_dia_only <- c("1004564", "1010661", "1010781", "1019964", "1020076", "1021087", "1023338", "1023616", "1024589", "1026233", "1026662", "B100", "EF7", "FB14", "G44", "HD7")

c <- read.delim("~/hdall/results/clinical/clinical_data.tsv", na.strings=c("", "NA", "n/a", "n/d", " "))
c <- c[!(c$patient_id %in% c("E", "RN14046", "1025678", "1021186")),]

#----
# CLEAN UP DATA
#----

c$age_dia <- as.numeric(c$age_dia)
c$sex[c$sex!="m" & c$sex!="f"] <- NA
c$sex <- as.factor(as.character(c$sex))

cr <- c[c$cohort=="relapsing",]
cr$endpoint[cr$endpoint==""] <- NA
cr$endpoint <- as.factor(as.character(cr$endpoint))
cr$blasts_dia <- as.numeric(as.character(cr$blasts_dia))
cr$blasts_rel <- as.numeric(as.character(cr$blasts_rel))
cr$first_rem_months <- as.numeric(as.character(cr$first_rem_months))
cr$second_rem_months <- as.numeric(as.character(cr$second_rem_months))

m <- read.delim("~/hdall/results/reseq/filtered-variants.reseq.cosmic.normaf.tsv")
m <- m[m$status!="REJECT" & m$non_silent==T & m$freq_leu >= 0.1,]

m.exome <- read.delim("~/hdall/results/filtered-variants.cosmic.normaf.tsv")
m.exome <- m.exome[m.exome$status!="REJECT" & m.exome$freq_leu >= 0.1,]
m.exome.ns <- m.exome[m.exome$non_silent==T,]

#----
# ADD ATTRIBUTES TO CLINICAL TABLE
#----
cr$age_rel <- cr$age_dia+cr$first_rem_months/12

cr$crebbp.dia <- ifelse(as.character(cr$patient_id) %in% patients_rel_only, NA, as.character(cr$patient_id) %in% as.character(m[m$gene=="CREBBP" & m$sample=="rem_dia", "patient"]))
cr$crebbp.rel <- ifelse(as.character(cr$patient_id) %in% patients_dia_only, NA, as.character(cr$patient_id) %in% as.character(m[m$gene=="CREBBP" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]))
cr$crebbp.all <- as.character(cr$patient_id) %in% as.character(m[m$gene=="CREBBP", "patient"])

cr$kras.dia <- ifelse(as.character(cr$patient_id) %in% patients_rel_only, NA, as.character(cr$patient_id) %in% as.character(m[m$gene=="KRAS" & m$sample=="rem_dia", "patient"]))
cr$kras.rel <- ifelse(as.character(cr$patient_id) %in% patients_dia_only, NA, as.character(cr$patient_id) %in% as.character(m[m$gene=="KRAS" & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]))
cr$kras.all <- as.character(cr$patient_id) %in% as.character(m[m$gene=="KRAS", "patient"])

cr$ras.dia <- ifelse(as.character(cr$patient_id) %in% patients_rel_only, NA, as.character(cr$patient_id) %in% as.character(m[(m$gene=="KRAS" | m$gene=="NRAS" | m$gene=="PTPN11" | m$gene=="FLT3") & m$sample=="rem_dia", "patient"]))
cr$ras.rel <- ifelse(as.character(cr$patient_id) %in% patients_dia_only, NA, as.character(cr$patient_id) %in% as.character(m[(m$gene=="KRAS" | m$gene=="NRAS" | m$gene=="PTPN11" | m$gene=="FLT3") & (m$sample=="rem_rel" | m$sample=="rem_rel2"), "patient"]))
cr$ras.all <- as.character(cr$patient_id) %in% as.character(m[m$gene=="KRAS" | m$gene=="NRAS" | m$gene=="PTPN11" | m$gene=="FLT3", "patient"])

cr$ras.crebbp.rel <- cr$crebbp.rel & cr$ras.rel

cr$death <- as.factor(cr$endpoint=="2nd rel + death" | cr$endpoint=="death")
cr$sec_rel <- as.factor(cr$endpoint=="2nd rel" | cr$endpoint=="2nd rel + death")

cr$num.mut.dia <- ifelse(as.character(cr$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome$patient==as.character(x) & m.exome$sample=="rem_dia") }), NA)
cr$num.mut.rel <- ifelse(as.character(cr$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome$patient==as.character(x) & m.exome$sample=="rem_rel") }), NA)
cr$num.mut.dia.ns <- ifelse(as.character(cr$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome.ns$patient==as.character(x) & m.exome.ns$sample=="rem_dia") }), NA)
cr$num.mut.rel.ns <- ifelse(as.character(cr$patient_id) %in% patients_exome, sapply(c$patient_id, function(x) { sum(m.exome.ns$patient==as.character(x) & m.exome.ns$sample=="rem_rel") }), NA)

#boxplot(first_rem_months ~ sex, data=cr, na.action=na.exclude, outline=F, names=c("female", "male"))
#stripchart(first_rem_months ~ sex, data=cr, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=c("red", "blue"), add=T)

pdf("~/hdall/results/clinical/clinical.pdf")

#par(mfrow=c(3,3), mar=c(2.5,4,1.5,2))

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