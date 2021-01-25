
## Code for generating Figure 3 from the main manuscript.
library(DelayedSurvFit)
data("NivoIpili")
head(NivoIpili)


## First look at median survival for OS in the two treatment groups.
## Should be 17.3 vs. 15.0 months
nivo.surv <- survfit(Surv(times, evnt) ~ trt, data = NivoIpili)
print(nivo.surv)

### Now, fit a single-crossing constrained estimate of the arm-specific survival curves
ds.obj <- DelayedSurvFit(times=NivoIpili$times, events=NivoIpili$evnt, 
                         trt=NivoIpili$trt, inner.iter=100)

#postscript("~/OSPD1_survival.eps", horizontal=FALSE, width=6.5, height=4.5)
par(mfrow=c(1,2), mar=c(3.9, 3.9, 1.5, .5))
plot(ds.obj$times, ds.obj$surv0, type="n", las=1, xlab="Months", ylab="Survival Probability",
     xlim=c(0, max(ds.obj$times)), ylim=c(0,1), cex.axis=0.8, cex.lab=0.9)
lines(c(0,ds.obj$times), c(1, ds.obj$surv0), type="s", lwd=2)
lines(c(0,ds.obj$times), c(1, ds.obj$surv1), type="s", col="red", lwd=2)
abline(v=ds.obj$theta, lwd=2)
legend("topright", legend=c("Nivolumab + ipilimumab", "Chemotherapy"), col=c("red", "black"), lwd=3,
       bty='n', cex=0.7)
text(x=26, y=.82, labels=expression(paste(hat(theta)[sc], " = ", 7.36)), cex=0.92)
arrows(19.25, .8, 7.36, 0.72, lwd=2, length=0.15)
text(x=3, y=.1, labels=expression(paste(hat(gamma)[sc], " = ", 1)), cex=0.8)

SS0 <- SurvFn(ds.obj, arm=0)
SS1 <- SurvFn(ds.obj, arm=1)
Stheta1 <- SS1(ds.obj$theta + 1e-4)
Stheta0 <- SS0(ds.obj$theta + 1e-4)
ind <- ds.obj$times >= ds.obj$theta

plot(ds.obj$times[ind], ds.obj$surv0[ind]/Stheta0, type="n", las=1, xlab="Months", ylab="Survival Prob.",
     xlim=c(ds.obj$theta, max(ds.obj$times)), ylim=c(0,1), cex.axis=0.8, cex.lab=0.9)
lines(ds.obj$times[ind], ds.obj$surv0[ind]/Stheta0, type="s", lwd=2)
lines(ds.obj$times[ind], ds.obj$surv1[ind]/Stheta1, type="s", col="red", lwd=2)
legend("topright", legend=c("Nivolumab + ipilimumab", "Chemotherapy"), col=c("red", "black"), lwd=3,
       bty='n', cex=0.7)
#dev.off()


## Conditional survival probabilities mentioned in paper
round(SS1(24)/SS1(ds.obj$theta), 2)
round(SS0(24)/SS0(ds.obj$theta), 2)
round(SS1(36)/SS1(ds.obj$theta), 2)
round(SS0(36)/SS0(ds.obj$theta), 2)




