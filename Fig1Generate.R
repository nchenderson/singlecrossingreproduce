
######################################################
## Code used to generate Figure 1 from the main manuscript.
## This figure shows the four possible survival profiles
######################################################

tt <- seq(0, 6, length.out=1000)
Sweibull <- function(x, alpha, lambda) {
    neglogS <- (x/lambda)^alpha 
    return(exp(-neglogS))
}

SS1a <- Sweibull(tt, alpha=1.5, lambda=6)
SS0a <- Sweibull(tt, alpha=1.5, lambda=4)
SS1b <- Sweibull(tt, alpha=2, lambda=4)
SS0b <- Sweibull(tt, alpha=.8, lambda=7)

ytxt <- -.2
ff <- function(x) {
   ans <- Sweibull(x, alpha=2, lambda=4) - Sweibull(x, alpha=.8, lambda=7)
   return(ans)
}
uniroot(ff, interval=c(.1, 6))


## Generate Figures: one black and white, one color.

#postscript("~/FourPossibleProfiles_color.eps", horizontal=FALSE, width=6, height=4.5)
par(mfrow=c(2,2), mar=c(2.6, 3.9, .5, .5))
plot(tt, SS1a, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n',
     cex.axis=0.9)
lines(tt, SS1a, col="red", lwd=2)
lines(tt, SS0a, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), col=c("red", "black"), bty='n', lwd=2)
legend("topright", legend=c(expression(paste(gamma," = 1")), expression(paste(theta," = 0"))), bty="n", cex=1.1)
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=.9)

plot(tt, SS1b, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1b, lwd=2)
lines(tt, SS0b, col="red", lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), col=c("red", "black"), bty='n', lwd=2)
legend("topright", legend=c(expression(paste(gamma," = 1")), expression(paste(theta," = 2.75"))), bty="n", cex=1.1)
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)

plot(tt, SS1a, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1a, lwd=2)
lines(tt, SS0a, lwd=2, col="red")
legend("bottomleft", legend=c("Treatment", "Control"), col=c("red", "black"), bty='n', lwd=2)
legend("topright", legend=c(expression(paste(gamma," = -1")), expression(paste(theta," = 0"))), bty="n", cex=1.1)
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)

plot(tt, SS1b, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1b, col="red", lwd=2)
lines(tt, SS0b, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), col=c("red", "black"), bty='n', lwd=2)
legend("topright", legend=c(expression(paste(gamma," = -1")), expression(paste(theta," = 2.75"))), bty="n", cex=1.1)
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)
#dev.off()


#postscript("~/FourPossibleProfiles.eps", horizontal=FALSE, width=6, height=4.5)
par(mfrow=c(2,2), mar=c(2.6, 3.9, .5, .5))
plot(tt, SS1a, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n',
     cex.axis=0.9)
lines(tt, SS1a, lwd=2)
lines(tt, SS0a, lty=2, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), lty=c(1,2), bty='n', lwd=1)
legend("topright", legend=expression(paste(gamma," = 1")), bty="n")
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=.9)

plot(tt, SS1b, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1b, lwd=2)
lines(tt, SS0b, lty=2, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), lty=c(1,2), bty='n', lwd=1)
legend("topright", legend=expression(paste(gamma," = -1")), bty="n")
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)

plot(tt, SS1a, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1a, lwd=2, lty=2)
lines(tt, SS0a, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), lty=c(1,2), bty='n', lwd=1)
legend("topright", legend=expression(paste(gamma," = -1")), bty="n")
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)

plot(tt, SS1b, type="n", ylim=c(0, 1), las=1, xlab="Time", ylab="Survival Probability", xaxt='n', cex.axis=0.9)
lines(tt, SS1b, lty=2, lwd=2)
lines(tt, SS0b, lwd=2)
legend("bottomleft", legend=c("Treatment", "Control"), lty=c(1,2), bty='n', lwd=1)
legend("topright", legend=expression(paste(gamma," = 1")), bty="n")
axis(1, at=0:6, labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=0.8)
text(x=3, y=ytxt, labels="Time", xpd=NA, cex=0.9)
#dev.off()

