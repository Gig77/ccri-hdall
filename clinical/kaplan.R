options(warn=1)
library(survival)
library(cmprsk)

rm(list=ls())
c <- read.delim("~/hdall/results/clinical/clinical_data.processed.tsv")

# ==========================================================================================================

pdf("~/hdall/results/clinical/kaplan-figure-paper.pdf", width=10, height=13.3)

par(mfrow=c(4,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

# --- MRD ---

form <- Surv(time=second_rem_months, dead)~mrd_risk_rel
fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 1), xlab="months", ylab="pOS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
#fit.overall <- survfit(Surv(time=c$second_rem_months[!is.na(c$mrd_risk_rel)], c$dead[!is.na(c$mrd_risk_rel)])~1)
#lines(fit.overall, col="black", lty=1, conf.int=F)
#text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="black", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "black", "red"), lty=c(1, 1, 1), box.lwd=0.5)
legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "red"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~mrd_risk_rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 1), xlab="months", ylab="pEFS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
#fit.overall <- survfit(Surv(time=c$second_rem_months[!is.na(c$mrd_risk_rel)], c$had_second_event_after_first[!is.na(c$mrd_risk_rel)])~1)
#lines(fit.overall, col="black", lty=1, conf.int=F)
#text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="black", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "black", "red"), lty=c(1, 1, 1), box.lwd=0.5)
legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "red"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

#RFS; do not use because of competitive risk (death); use cumulative incidence instead
#rm(form); form <- Surv(time=second_rem_months, sec_rel)~mrd_risk_rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("red", "blue"), lty=c(1, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
#axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
#text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
#fit.overall <- survfit(Surv(time=c$second_rem_months[!is.na(c$mrd_risk_rel)], c$sec_rel[!is.na(c$mrd_risk_rel)])~1) 
#lines(fit.overall, col="black", lty=1, conf.int=F)
#text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="black", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "black", "red"), lty=c(1, 1, 1), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

#cumulative incidence
fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$mrd_risk_rel, cencode="none")
plot(fit[1:2], xlab="months", ylab="CIR", curvlab=c("MRD HR", "MRD SR"), col=c("red", "blue"), wh=c(-100,-100), yaxt='n', lty=c(1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["HR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["HR 2nd rel",], timepoints(fit, 60)$var["HR 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["SR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["SR 2nd rel",], timepoints(fit, 60)$var["SR 2nd rel",]), col="blue", adj=0, cex=1.3)
#fit.overall <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, cencode="none") 
#tp <- c(0, unique(c$second_rem_months[!is.na(c$second_rem_months) & !is.na(c$mrd_risk_rel) & c$second_event_after_first_relapse=="2nd rel"]), max(c[!is.na(c$second_rem_months) & !is.na(c$mrd_risk_rel) & !is.na(c$second_event_after_first_relapse), "second_rem_months"]))
#tp <- tp[order(tp)]
#lines(tp, timepoints(fit.overall, tp)$est[1,], col="black", lty=1, type="s")
#text(60, timepoints(fit.overall, 60)$est[1,]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit.overall, 60)$est[1,], timepoints(fit.overall, 60)$var[1,]), col="black", adj=0, cex=1.3)
box()
#legend("topright", c(sprintf("MRD HR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="HR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="HR",c("second_rem_months", "second_event_after_first_relapse")]))),
#				    sprintf("overall (%d/%d)", sum(complete.cases(c[!is.na(c$mrd_risk_rel) & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[!is.na(c$mrd_risk_rel),c("second_rem_months", "second_event_after_first_relapse")]))),
#				    sprintf("MRD SR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="SR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="SR",c("second_rem_months", "second_event_after_first_relapse")])))), 
#			lwd=c(1,1), col=c("red", "black", "blue"), lty=c(1, 1, 1), box.lwd=0.5)
legend("topright", c(sprintf("MRD HR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="HR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="HR",c("second_rem_months", "second_event_after_first_relapse")]))),
					sprintf("MRD SR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="SR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="SR",c("second_rem_months", "second_event_after_first_relapse")])))), 
			lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

# --- RAS PW ---

rm(form); form <- Surv(time=second_rem_months, dead)~ras.relapsing.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pOS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.relapsing.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pEFS", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

#RFS; do not use because of competitive risk (death); use cumulative incidence instead
#rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n', cex.axis=1.3, cex.lab=1.5)
#axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
#text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0, cex=1.3)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$ras.relapsing.rel, cencode="none")
plot(fit[1:2], xlab="months", ylab="CIR", curvlab=c("wtRAS", "mRAS"), col=c("blue", "red"), wh=c(-100,-100), yaxt='n', lty=c(1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["TRUE 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["TRUE 2nd rel",], timepoints(fit, 60)$var["TRUE 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["FALSE 2nd rel",]-0.05, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["FALSE 2nd rel",], timepoints(fit, 60)$var["FALSE 2nd rel",]), col="blue", adj=0, cex=1.3)
box()
legend("topright", c(sprintf("mRAS (%d/%d)", sum(complete.cases(c[c$ras.relapsing.rel & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$ras.relapsing.rel, c("second_rem_months", "second_event_after_first_relapse")]))), 
				   sprintf("wtRAS (%d/%d)", sum(complete.cases(c[!c$ras.relapsing.rel & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[!c$ras.relapsing.rel, c("second_rem_months", "second_event_after_first_relapse")])))), 
		lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

# --- CREBBP ---

#rm(form); form <- Surv(time=second_rem_months, dead)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pOS", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0)

#rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pEFS", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0)

#rm(form); form <- Surv(time=second_rem_months, sec_rel)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(1, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(62, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 1), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0)

# --- RAS PW BY GENE ---

c$ras.pw.gene.rel <- as.character(c$ras.pw.gene.rel)
#c$ras.pw.gene.rel[c$ras.pw.gene.rel=="FLT3"] <- "WT"  # lump FLT3 into WT group for now (only few cases)

rm(form); form <- Surv(time=second_rem_months, dead)~ras.pw.gene.rel
rm(fit); fit <- survfit(form, data=c[c$ras.pw.gene.rel!="FLT3",]) 
plot(fit, col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), xlab="months", ylab="pOS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="#117733", adj=0, cex=1.3)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(1, 1, 1, 1), col=c("#117733", "blue", "black", "red"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$chisq, length(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.pw.gene.rel
rm(fit); fit <- survfit(form, data=c[c$ras.pw.gene.rel!="FLT3",]) 
plot(fit, col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), xlab="months", ylab="pEFS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="#117733", adj=0, cex=1.3)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(1, 1, 1, 1), col=c("#117733", "blue", "black", "red"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$chisq, length(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$n) - 1))))), adj=0, cex=1.3)

#RFS; do not use because of competitive risk (death); use cumulative incidence instead
#rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.pw.gene.rel
#rm(fit); fit <- survfit(form, data=c[c$ras.pw.gene.rel!="FLT3",]) 
#plot(fit, col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
#axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
#text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
#text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="#117733", adj=0, cex=1.3)
#text(60, fit[4]$surv[max(which(fit[4]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(1, 1, 1, 1), col=c("#117733", "blue", "black", "red"), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$chisq, length(survdiff(form, data=c[c$ras.pw.gene.rel!="FLT3",])$n) - 1))))), adj=0, cex=1.3)

fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$ras.pw.gene.rel, cencode="none", subset=c$ras.pw.gene.rel!="FLT3")
plot(fit[1:4], xlab="months", ylab="CIR", curvlab=c("mKRAS", "mNRAS", "mPTPN11", "wtRAS"), col=c("red", "blue", "#117733", "black"), wh=c(-100,-100), yaxt='n', lty=c(1, 1, 1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["KRAS 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["KRAS 2nd rel",], timepoints(fit, 60)$var["KRAS 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["NRAS 2nd rel",]-0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["NRAS 2nd rel",], timepoints(fit, 60)$var["NRAS 2nd rel",]), col="blue", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["PTPN11 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["PTPN11 2nd rel",], timepoints(fit, 60)$var["PTPN11 2nd rel",]), col="#117733", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["WT 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["WT 2nd rel",], timepoints(fit, 60)$var["WT 2nd rel",]), col="black", adj=0, cex=1.3)
box()
legend("topright", c(sprintf("mKRAS (%d/%d)", sum(complete.cases(c[c$ras.pw.gene.rel=="KRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$ras.pw.gene.rel=="KRAS", c("second_rem_months", "second_event_after_first_relapse")]))),
				    sprintf("wtRAS (%d/%d)", sum(complete.cases(c[c$ras.pw.gene.rel=="WT" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$ras.pw.gene.rel=="WT", c("second_rem_months", "second_event_after_first_relapse")]))),
				    sprintf("mNRAS (%d/%d)", sum(complete.cases(c[c$ras.pw.gene.rel=="NRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$ras.pw.gene.rel=="NRAS", c("second_rem_months", "second_event_after_first_relapse")]))), 
					sprintf("mPTPN11 (%d/%d)", sum(complete.cases(c[c$ras.pw.gene.rel=="PTPN11" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$ras.pw.gene.rel=="PTPN11", c("second_rem_months", "second_event_after_first_relapse")]))) 
					), 
					lwd=c(1,1), col=c("red", "black", "blue", "#117733"), lty=c(1, 1, 1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

# --- RAS/CREBBP co-occurrence ---

c$kras.crebbp.rel <- NA
c$kras.crebbp.rel[c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP and mKRAS"
c$kras.crebbp.rel[c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "mCREBBP or mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP or mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP and wtKRAS"

rm(form); form <- Surv(time=second_rem_months, dead)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "black"), lty=c(1, 1, 1), xlab="months", ylab="pOS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="black", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mCREBBP and mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("mCREBBP or mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), 
				sprintf("wtCREBBP and wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n)), 
		lty=c(1, 1, 1), col=c("red", "blue", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "black"), lty=c(1, 1, 1), xlab="months", ylab="pEFS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="black", adj=0, cex=1.3)
legend("bottomright", c(sprintf("mCREBBP and mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("mCREBBP or mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), 
				sprintf("wtCREBBP and wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n)), 
		lty=c(1, 1, 1), col=c("red", "blue", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$kras.crebbp.rel, cencode="none")
plot(fit[1:3], xlab="months", ylab="CIR", curvlab=c("mCREBBP and mKRAS", "mCREBBP or mKRAS", "wtCREBBP and wtKRAS"), col=c("red", "blue", "black"), wh=c(-100,-100), yaxt='n', lty=c(1, 1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["mCREBBP and mKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["mCREBBP and mKRAS 2nd rel",], timepoints(fit, 60)$var["mCREBBP and mKRAS 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["mCREBBP or mKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["mCREBBP or mKRAS 2nd rel",], timepoints(fit, 60)$var["mCREBBP or mKRAS 2nd rel",]), col="blue", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["wtCREBBP and wtKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["wtCREBBP and wtKRAS 2nd rel",], timepoints(fit, 60)$var["wtCREBBP and wtKRAS 2nd rel",]), col="black", adj=0, cex=1.3)
box()
legend("topright", c(sprintf("mCREBBP and mKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP and mKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP and mKRAS", c("second_rem_months", "second_event_after_first_relapse")]))), 
				sprintf("mCREBBP or mKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP or mKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP or mKRAS", c("second_rem_months", "second_event_after_first_relapse")]))),
				sprintf("wtCREBBP and wtKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP and wtKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP and wtKRAS", c("second_rem_months", "second_event_after_first_relapse")]))) 
		), 
		lwd=c(1,1), col=c("red", "blue", "black"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

c$mkras.wtcrebbp.relapsing.rel <- c$kras.relapsing.rel & !c$crebbp.relapsing.rel
c$wtkras.mcrebbp.relapsing.rel <- !c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$mkras.mcrebbp.relapsing.rel <- c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$wtkras.wtcrebbp.relapsing.rel <- !c$kras.relapsing.rel & !c$crebbp.relapsing.rel

dev.off()

# ==========================================================================================================

pdf("~/hdall/results/clinical/kaplan-ras-gene-vs-ras-wt-pairwise.pdf", width=10, height=13.3)

par(mfrow=c(4,3),mar=c(4,4,1,1))

# --- KRAS vs. WT ---

c$ras.pw.gene.rel.pairwise <- ifelse(c$ras.pw.gene.rel %in% c("NRAS", "PTPN11", "FLT3"), NA, c$ras.pw.gene.rel)

rm(form); form <- Surv(time=second_rem_months, dead)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

# --- NRAS vs. WT ---

c$ras.pw.gene.rel.pairwise <- ifelse(c$ras.pw.gene.rel %in% c("KRAS", "PTPN11", "FLT3"), NA, c$ras.pw.gene.rel)

rm(form); form <- Surv(time=second_rem_months, dead)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mNRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mNRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mNRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

# --- PTPN11 vs. WT ---

c$ras.pw.gene.rel.pairwise <- ifelse(c$ras.pw.gene.rel %in% c("KRAS", "NRAS", "FLT3"), NA, c$ras.pw.gene.rel)

rm(form); form <- Surv(time=second_rem_months, dead)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mPTPN11 (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mPTPN11 (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.pw.gene.rel.pairwise
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black"), lty=c(1, 2), xlab="months", ylab="pOS", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("wtRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("mPTPN11 (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(2, 1), col=c("black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

dev.off()

stop("OK")


# --- old stuff ----

pdf("~/hdall/results/clinical/kaplan-mutations-at-relapse.pdf")

plot(survfit(Surv(time=c$second_rem_months, c$dead)~1), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)

form <- Surv(time=second_rem_months, dead)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pRFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~nras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(80, 0.2, c("wtNRAS", "mNRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~nras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS", conf.int=F)
legend(80, 0.2, c("wtNRAS", "mNRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~ptpn11.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(80, 0.2, c("wtPTPN11", "mPTPN11"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ptpn11.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS", conf.int=F)
legend(80, 0.2, c("wtPTPN11", "mPTPN11"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)


rm(form); form <- Surv(time=second_rem_months, dead)~crebbp.and.kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(50, 0.2, c("wtCREBBP or wtKRAS", "mCREBBP and mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~crebbp.and.kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS", conf.int=F)
legend(50, 0.2, c("wtCREBBP or wtKRAS", "mCREBBP and mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~ikzf1.del.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(60, 0.2, c("wtIKZF1", "IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ikzf1.del.rel
plot(survfit(form, data=c), col=c("red", "blue"), xlab="months", ylab="pEFS", conf.int=F)
legend(80, 0.2, c("wtIKZF1", "IKZF1-del"), lwd=c(1,1), col=c("red", "blue"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS", conf.int=F)
legend(60, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS", conf.int=F)
legend(50, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pRFS after 1st relapse", conf.int=F)
legend(60, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

dev.off()

# --- RAS/CREBBP co-occurrence, split into all four combinations ---

pdf("~/hdall/results/clinical/kaplan-kras-crebbp.pdf", width=10, height=6.65)
par(mfrow=c(2,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

c$kras.crebbp.rel <- NA
c$kras.crebbp.rel[c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP+wtKRAS"
c$kras.crebbp.rel[c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP+mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+wtKRAS"

rm(form); form <- Surv(time=second_rem_months, dead)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), xlab="months", ylab="pOS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="#117733", adj=0, cex=1.3)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0, cex=1.3)
legend(40, 0.5, c(sprintf("mCREBBP+mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("mCREBBP+wtKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n),
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), 
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n)), 
		lty=c(1, 1, 1, 1), col=c("red", "blue", "#117733", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), xlab="months", ylab="pEFS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="#117733", adj=0, cex=1.3)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0, cex=1.3)
legend(40, 0.5, c(sprintf("mCREBBP+mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("mCREBBP+wtKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n),
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), 
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n)), 
		lty=c(1, 1, 1, 1), col=c("red", "blue", "#117733", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$kras.crebbp.rel, cencode="none")
plot(fit[1:4], xlab="months", ylab="CIR", curvlab=c("mCREBBP+mKRAS", "mCREBBP+wtKRAS", "wtCREBBP+mKRAS", "wtCREBBP+wtKRAS"), col=c("red", "blue", "#117733", "black"), wh=c(-100,-100), yaxt='n', lty=c(1, 1, 1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["mCREBBP+mKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["mCREBBP+mKRAS 2nd rel",], timepoints(fit, 60)$var["mCREBBP+mKRAS 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["mCREBBP+wtKRAS 2nd rel",]-0.05, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["mCREBBP+wtKRAS 2nd rel",], timepoints(fit, 60)$var["mCREBBP+wtKRAS 2nd rel",]), col="blue", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["wtCREBBP+mKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["wtCREBBP+mKRAS 2nd rel",], timepoints(fit, 60)$var["wtCREBBP+mKRAS 2nd rel",]), col="#117733", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["wtCREBBP+wtKRAS 2nd rel",]-0.05, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["wtCREBBP+wtKRAS 2nd rel",], timepoints(fit, 60)$var["wtCREBBP+wtKRAS 2nd rel",]), col="black", adj=0, cex=1.3)
box()
legend("topright", c(sprintf("mCREBBP+mKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP+mKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP+mKRAS", c("second_rem_months", "second_event_after_first_relapse")]))),
				sprintf("mCREBBP+wtKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP+wtKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP+wtKRAS", c("second_rem_months", "second_event_after_first_relapse")]))), 
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+mKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+mKRAS", c("second_rem_months", "second_event_after_first_relapse")]))), 
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+wtKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+wtKRAS", c("second_rem_months", "second_event_after_first_relapse")])))), 
		lwd=c(1,1), col=c("red", "blue", "#117733", "black"), lty=c(1, 1, 1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)


#---


c$kras.crebbp.rel <- NA
c$kras.crebbp.rel[c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP"
c$kras.crebbp.rel[c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP"
c$kras.crebbp.rel[!c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+wtKRAS"

rm(form); form <- Surv(time=second_rem_months, dead)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black", "blue"), lty=c(1, 1, 1), xlab="months", ylab="pOS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0, cex=1.3)
legend("right", c(sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), 
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(1, 1, 1), col=c("red", "blue", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "black", "blue"), lty=c(1, 1, 1), xlab="months", ylab="pEFS", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0, cex=1.3)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="black", adj=0, cex=1.3)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0, cex=1.3)
legend("right", c(sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n),
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n),  
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(1, 1, 1), col=c("red", "blue", "black"), box.lwd=0.5)
text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

#RFS; do not use because of competitive risk (death); use cumulative incidence instead
#rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.crebbp.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("black", "red", "blue"), lty=c(1, 1, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, yaxt='n', cex.axis=1.3, cex.lab=1.5)
#axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
#text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="black", adj=0, cex=1.3)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0, cex=1.3)
#text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.04, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0, cex=1.3)
#legend("bottomright", c(sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n), sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(1, 1, 1), col=c("blue", "black", "red"), box.lwd=0.5)
#text(1, 0, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1))))), adj=0, cex=1.3)

fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$kras.crebbp.rel, cencode="none")
plot(fit[1:3], xlab="months", ylab="CIR", curvlab=c("mCREBBP", "wtCREBBP+mKRAS", "wtCREBBP+wtKRAS"), col=c("red", "black", "blue"), wh=c(-100,-100), yaxt='n', lty=c(1, 1, 1), xlim=c(0,125), cex.axis=1.3, cex.lab=1.5)
axis(side=2, at=seq(0, 1, by=0.1), cex.axis=1.3)
text(60, timepoints(fit, 60)$est["mCREBBP 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["mCREBBP 2nd rel",], timepoints(fit, 60)$var["mCREBBP 2nd rel",]), col="red", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["wtCREBBP+mKRAS 2nd rel",]-0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["wtCREBBP+mKRAS 2nd rel",], timepoints(fit, 60)$var["wtCREBBP+mKRAS 2nd rel",]), col="black", adj=0, cex=1.3)
text(60, timepoints(fit, 60)$est["wtCREBBP+wtKRAS 2nd rel",]+0.04, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["wtCREBBP+wtKRAS 2nd rel",], timepoints(fit, 60)$var["wtCREBBP+wtKRAS 2nd rel",]), col="blue", adj=0, cex=1.3)
box()
legend("topright", c(sprintf("mCREBBP (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="mCREBBP", c("second_rem_months", "second_event_after_first_relapse")]))), 
				sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+wtKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+wtKRAS", c("second_rem_months", "second_event_after_first_relapse")]))),
				sprintf("wtCREBBP+mKRAS (%d/%d)", sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+mKRAS" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$kras.crebbp.rel=="wtCREBBP+mKRAS", c("second_rem_months", "second_event_after_first_relapse")]))) 
		), 
		lwd=c(1,1), col=c("red", "blue", "black"), lty=c(1, 1, 1), box.lwd=0.5)
text(1, 1, substitute(italic(P) == p, list(p=gsub("0\\.", "\\.", sprintf("%.2g", fit$Tests["2nd rel", "pv"])))), adj=0, cex=1.3)

dev.off()
