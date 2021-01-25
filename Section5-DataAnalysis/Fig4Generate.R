
## Code for generating Figure 4 from the main manuscript.
library(DelayedSurvFit)
data("NivoIpili")
head(NivoIpili)

ds.obj <- DelayedHazFit(times=NivoIpili$times, events=NivoIpili$evnt, 
                        trt=NivoIpili$trt, inner.iter=100)

#postscript("~/OSPD1_hazard.eps", horizontal=FALSE, width=6.5, height=4.5)
par(mfrow=c(1,2), mar=c(3.9, 3.9, 1.5, .5))
plot(ds.obj$times, ds.obj$haz0, type="n", las=1, xlab="Months", ylab="hazard",
     xlim=c(0, max(ds.obj$times)), ylim=c(0,.072), cex.axis=0.8, cex.lab=0.9)
ntimes <- length(ds.obj$times)
for(k in 1:ntimes) {
  if(ds.obj$haz0[k] > ds.obj$haz1[k]) {
     lines(c(ds.obj$times[k], ds.obj$times[k]), c(0, ds.obj$haz0[k]), lwd=2)
     lines(c(ds.obj$times[k], ds.obj$times[k]), c(0, ds.obj$haz1[k]), col="red", lwd=2)
  } else {
     lines(c(ds.obj$times[k], ds.obj$times[k]), c(0, ds.obj$haz1[k]), col="red", lwd=2)
     lines(c(ds.obj$times[k], ds.obj$times[k]), c(0, ds.obj$haz0[k]), lwd=2)
  }
}
abline(v=ds.obj$theta, lwd=2, lty=2)
legend("top", legend=c("Nivolumab + ipilimumab","Chemotherapy"), col=c("red", "black"), lwd=3,
       bty='n', cex=0.7)
text(x=16.5, y=.055, labels=expression(paste(hat(theta)[sc], " = ", 2.4)))
arrows(10.0, .054, 2.4, 0.045, lwd=2, length=0.15)

tmp_fit0 <- lowess(ds.obj$times, ds.obj$haz0, f=0.33)
tmp_fit1 <- lowess(ds.obj$times, ds.obj$haz1, f=0.33)
plot(tmp_fit0$x, tmp_fit0$y, type="n", las=1, xlab="Months", ylab="hazard",
     xlim=c(0, 36), ylim=c(0,.03), cex.axis=0.8, cex.lab=0.9)
lines(tmp_fit0$x, tmp_fit0$y, lwd=2)
lines(tmp_fit1$x, tmp_fit1$y, col="red", lwd=2)
legend("topright", legend=c("Nivolumab + ipilimumab","Chemotherapy"), col=c("red", "black"), lwd=3,
       bty='n', cex=0.7)
#dev.off()


## pre- and post-crossing hazard ratios referred to in the main manuscript
tmp <- AvgHazardRatios(ds.obj)
round(tmp$precross.avg, 2)
round(tmp$postcross.avg, 2)











