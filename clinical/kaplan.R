options(warn=1)
library(survival)
library(cmprsk)

c <- read.delim("~/hdall/results/clinical/clinical_data.processed.tsv")

pdf("~/hdall/results/clinical/kaplan-figure-paper.pdf", width=10, height=13.3)

par(mfrow=c(4,3),mar=c(4,4,1,1))

# --- MRD ---

rm(form); form <- Surv(time=second_rem_months, dead)~mrd_risk_rel
fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 5), xlab="months", ylab="pOS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
fit.overall <- survfit(Surv(time=c$second_rem_months, c$dead)~1)
lines(fit.overall, col="darkgray", lty=2, conf.int=F)
text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="darkgray", adj=0)
legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "darkgray", "red"), lty=c(5, 2, 1))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~mrd_risk_rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 5), xlab="months", ylab="pEFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
fit.overall <- survfit(Surv(time=c$second_rem_months, c$had_second_event_after_first)~1)
lines(fit.overall, col="darkgray", lty=2, conf.int=F)
text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="darkgray", adj=0)
legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "darkgray", "red"), lty=c(5, 2, 1))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

#RFS; do not use because of competitive risk (death); use cumulative incidence instead
rm(form); form <- Surv(time=second_rem_months, sec_rel)~mrd_risk_rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue"), lty=c(1, 5), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
fit.overall <- survfit(Surv(time=c$second_rem_months, c$sec_rel)~1) 
lines(fit.overall, col="darkgray", lty=2, conf.int=F)
text(60, fit.overall$surv[max(which(fit.overall$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit.overall$surv[max(which(fit.overall$time<=60))], fit.overall$std.err[max(which(fit.overall$time<=60))]), col="darkgray", adj=0)
legend("bottomright", c(sprintf("MRD SR (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("overall (%d/%d)", sum(fit.overall$n.event), fit.overall$n), sprintf("MRD HR (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("blue", "darkgray", "red"), lty=c(5, 2, 1))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

#cumulative incidence
#fit <- cuminc(c$second_rem_months, c$second_event_after_first_relapse, c$mrd_risk_rel, cencode="none")
#plot(fit[1:2], xlab="months", ylab="cumulative incidence 2nd relapse", curvlab=c("MRD HR", "MRD SR"), col=c("red", "blue"), wh=c(-100,-100), yaxt='n', lty=c(1, 2), xlim=c(0,125))
#axis(side=2, at=seq(0, 1, by=0.1))
#text(60, timepoints(fit, 60)$est["HR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["HR 2nd rel",], timepoints(fit, 60)$var["HR 2nd rel",]), col="red", adj=0)
#text(60, timepoints(fit, 60)$est["SR 2nd rel",]+0.03, sprintf("%.2f, SE=%.2f", timepoints(fit, 60)$est["SR 2nd rel",], timepoints(fit, 60)$var["SR 2nd rel",]), col="blue", adj=0)
#legend("topleft", c(sprintf("MRD HR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="HR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="HR",c("second_rem_months", "second_event_after_first_relapse")]))), 
#				    sprintf("MRD SR (%d/%d)", sum(complete.cases(c[c$mrd_risk_rel=="SR" & c$second_event_after_first_relapse=="2nd rel", "second_rem_months"])),  sum(complete.cases(c[c$mrd_risk_rel=="SR",c("second_rem_months", "second_event_after_first_relapse")])))), 
#			lwd=c(1,1), col=c("red", "blue"), lty=c(1, 2))
#text(120, 1, sprintf("p=%.2g", fit$Tests["2nd rel", "pv"]), adj=1)

# --- RAS PW ---

rm(form); form <- Surv(time=second_rem_months, dead)~ras.relapsing.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pOS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.relapsing.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pEFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.relapsing.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
legend("bottomright", c(sprintf("mRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

# --- CREBBP ---

#rm(form); form <- Surv(time=second_rem_months, dead)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pOS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
#text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

#rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pEFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
#text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

#rm(form); form <- Surv(time=second_rem_months, sec_rel)~crebbp.relapsing.rel
#rm(fit); fit <- survfit(form, data=c) 
#plot(fit, col=c("blue", "red"), lty=c(5, 1), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, xlim=c(0,125), yaxt='n')
#axis(side=2, at=seq(0, 1, by=0.1))
#text(63, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="blue", adj=0)
#text(62, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
#legend("bottomright", c(sprintf("mCREBBP (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lwd=c(1,1), col=c("red", "blue"), lty=c(1, 5))
#text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

# --- RAS PW BY GENE ---

c$ras.pw.gene.rel <- as.character(c$ras.pw.gene.rel)
c$ras.pw.gene.rel[c$ras.pw.gene.rel=="FLT3"] <- "WT"  # lump FLT3 into WT group for now (only few cases)

rm(form); form <- Surv(time=second_rem_months, dead)~ras.pw.gene.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "orange", "black"), lty=c(1, 3, 4, 2), xlab="months", ylab="pOS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="orange", adj=0)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(4, 3, 2, 1), col=c("orange", "blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ras.pw.gene.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "orange", "black"), lty=c(1, 3, 4, 2), xlab="months", ylab="pEFS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="orange", adj=0)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(4, 3, 2, 1), col=c("orange", "blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~ras.pw.gene.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("red", "blue", "orange", "black"), lty=c(1, 3, 4, 2), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="red", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="blue", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]-0.05, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="orange", adj=0)
text(60, fit[4]$surv[max(which(fit[4]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[4]$surv[max(which(fit[4]$time<=60))], fit[4]$std.err[max(which(fit[4]$time<=60))]), col="black", adj=0)
legend("bottomright", c(sprintf("mPTPN11 (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mNRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n), sprintf("wtRAS (%d/%d)", sum(fit[4]$n.event), fit[4]$n), sprintf("mKRAS (%d/%d)", sum(fit[1]$n.event), fit[1]$n)), lty=c(4, 3, 2, 1), col=c("orange", "blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

# --- RAS/CREBBP co-occurrence ---

c$kras.crebbp.rel <- NA
c$kras.crebbp.rel[c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+mKRAS"
c$kras.crebbp.rel[!c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP"
c$kras.crebbp.rel[c$kras.relapsing.rel & c$crebbp.relapsing.rel] <- "mCREBBP"
c$kras.crebbp.rel[!c$kras.relapsing.rel & !c$crebbp.relapsing.rel] <- "wtCREBBP+wtKRAS"

rm(form); form <- Surv(time=second_rem_months, dead)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("black", "red", "blue"), lty=c(3, 1, 4), xlab="months", ylab="pOS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="black", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0)
legend("bottomright", c(sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n), sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(4, 3, 1), col=c("blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("black", "red", "blue"), lty=c(3, 1, 4), xlab="months", ylab="pEFS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="black", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0)
legend("bottomright", c(sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n), sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(4, 3, 1), col=c("blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.crebbp.rel
rm(fit); fit <- survfit(form, data=c) 
plot(fit, col=c("black", "red", "blue"), lty=c(3, 1, 4), xlab="months", ylab="pRFS after 1st relapse", conf.int=F, yaxt='n')
axis(side=2, at=seq(0, 1, by=0.1))
text(60, fit[1]$surv[max(which(fit[1]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[1]$surv[max(which(fit[1]$time<=60))], fit[1]$std.err[max(which(fit[1]$time<=60))]), col="black", adj=0)
text(60, fit[2]$surv[max(which(fit[2]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[2]$surv[max(which(fit[2]$time<=60))], fit[2]$std.err[max(which(fit[2]$time<=60))]), col="red", adj=0)
text(60, fit[3]$surv[max(which(fit[3]$time<=60))]+0.03, sprintf("%.2f, SE=%.2f", fit[3]$surv[max(which(fit[3]$time<=60))], fit[3]$std.err[max(which(fit[3]$time<=60))]), col="blue", adj=0)
legend("bottomright", c(sprintf("wtCREBBP+wtKRAS (%d/%d)", sum(fit[3]$n.event), fit[3]$n), sprintf("mCREBBP (%d/%d)", sum(fit[1]$n.event), fit[1]$n), sprintf("wtCREBBP+mKRAS (%d/%d)", sum(fit[2]$n.event), fit[2]$n)), lty=c(4, 3, 1), col=c("blue", "black", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

c$mkras.wtcrebbp.relapsing.rel <- c$kras.relapsing.rel & !c$crebbp.relapsing.rel
c$wtkras.mcrebbp.relapsing.rel <- !c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$mkras.mcrebbp.relapsing.rel <- c$kras.relapsing.rel & c$crebbp.relapsing.rel
c$wtkras.wtcrebbp.relapsing.rel <- !c$kras.relapsing.rel & !c$crebbp.relapsing.rel

dev.off()

stop("OK")


# --- old stuff ----

pdf("~/hdall/results/clinical/kaplan-mutations-at-relapse.pdf")

plot(survfit(Surv(time=c$second_rem_months, c$dead)~1), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)

form <- Surv(time=second_rem_months, dead)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pRFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtKRAS", "mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~nras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtNRAS", "mNRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~nras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtNRAS", "mNRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~ptpn11.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtPTPN11", "mPTPN11"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ptpn11.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtPTPN11", "mPTPN11"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)


rm(form); form <- Surv(time=second_rem_months, dead)~crebbp.and.kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(50, 0.2, c("wtCREBBP or wtKRAS", "mCREBBP and mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~crebbp.and.kras.relapsing.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(50, 0.2, c("wtCREBBP or wtKRAS", "mCREBBP and mKRAS"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~ikzf1.del.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(60, 0.2, c("wtIKZF1", "IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~ikzf1.del.rel
plot(survfit(form, data=c), col=c("red", "blue"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(80, 0.2, c("wtIKZF1", "IKZF1-del"), lwd=c(1,1), col=c("red", "blue"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, dead)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pOS after 1st relapse", conf.int=F)
legend(60, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, had_second_event_after_first)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pEFS after 1st relapse", conf.int=F)
legend(50, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

rm(form); form <- Surv(time=second_rem_months, sec_rel)~kras.or.ikzf1.rel
plot(survfit(form, data=c), col=c("blue", "red"), xlab="months", ylab="pRFS after 1st relapse", conf.int=F)
legend(60, 0.2, c("wtKRAS and wtIKZF1", "mKRAS or IKZF1-del"), lwd=c(1,1), col=c("blue", "red"))
text(1, 0, sprintf("p=%.2g", 1 - pchisq(survdiff(form, data=c)$chisq, length(survdiff(form, data=c)$n) - 1)), adj=0)

dev.off()
