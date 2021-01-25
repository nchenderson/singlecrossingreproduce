
### Code to generate Figure 2 from the main manuscript
### This figure plots the survival curves from each of the six simulation scenarios

Spwexp_point <- function(tt, rate=1, intervals, xmax) {
  n.int <- length(intervals)
  x <- y <- rep(0, n.int + 1)
  x[1] <- 0
  y[1] <- 0
  
  if(n.int > 1) {
    for(k in 2:n.int) {
      x[k] <- x[k-1] + intervals[k-1]
      y[k] <- y[k-1] + rate[k-1]*intervals[k-1]
    }
  }
  x[n.int + 1] <- xmax
  y[n.int + 1] <- y[n.int] + rate[n.int]*(xmax - sum(intervals[-n.int]))
  tmp <- approxfun(x, y)
  x.grid <- seq(0, xmax, length.out=1000)
  HH.x <- tmp(tt)
  return(exp(-HH.x))
}

Spwexp <- function(rate=1, intervals, xmax) {
  n.int <- length(intervals)
  x <- y <- rep(0, n.int + 1)
  x[1] <- 0
  y[1] <- 0
  
  if(n.int > 1) {
    for(k in 2:n.int) {
      x[k] <- x[k-1] + intervals[k-1]
      y[k] <- y[k-1] + rate[k-1]*intervals[k-1]
    }
  }
  x[n.int + 1] <- xmax
  y[n.int + 1] <- y[n.int] + rate[n.int]*(xmax - sum(intervals[-n.int]))
  tmp <- approxfun(x, y)
  x.grid <- seq(0, xmax, length.out=1000)
  HH <- tmp(x.grid)
  return(list(x=x.grid, y=exp(-HH), H=HH))
}



############################################################
# Simulation Study 1 (one survival curve dominates the other)

r01 <- c(.75, .3, .15, .1, .05)
r11 <- c(.5, .25, .125, .075, .025)
int1 <- c(rep(1, 4), Inf)

#############################################################
#### Second simulation study
## Have a very clear crossing

r12 <- c(.3, .25, .2, .1, .05, .05)
r02 <- c(.1, .15, .25, .2, .2, .2)
int2 <- c(rep(1, 5), Inf)

ff <- function(x) {
  ans <- Spwexp_point(x, rate=r12, intervals=int2, xmax=8) -  Spwexp_point(x, rate=r02, intervals=int2, xmax=8)
  return(ans)
}
cross2 <- uniroot(ff, interval=c(.1, 7.9))$root
cross2.y <- Spwexp_point(cross2, rate=r12, intervals=int2, xmax=8)

######################################################
### Third simulation study (Early Crossing)

r13 <- c(.35, .2, .15, .15, .15, .1, .1, .1, .1, .1, .1)
r03 <- c(.1, .25, .25, .25, .25, .2, .2, .2, .2, .2, .2)
int3 <- c(rep(1/2, 10), Inf)

ff <- function(x) {
  ans <- Spwexp_point(x, rate=r13, intervals=int3, xmax=8) -  Spwexp_point(x, rate=r03, intervals=int3, xmax=8)
  return(ans)
}
cross3 <- uniroot(ff, interval=c(.1, 7.9))$root
cross3.y <- Spwexp_point(cross3, rate=r13, intervals=int3, xmax=8)



######################################################
### Fourth simulation study (Early Crossing with almost no separation)

r14 <- c(.3, .2, .2, .15, .15, .1, .1, .1, .1, .1, .1)
r04 <- c(.275, .25, .25, .25, .25, .2, .2, .2, .2, .2, .2)
int4 <- c(rep(1/2, 10), Inf)

ff <- function(x) {
  ans <- Spwexp_point(x, rate=r14, intervals=int4, xmax=8) -  Spwexp_point(x, rate=r04, intervals=int4, xmax=8)
  return(ans)
}
cross4 <- uniroot(ff, interval=c(.1, 7.9))$root
cross4.y <- Spwexp_point(cross4, rate=r14, intervals=int4, xmax=8)

######################################################
### Fifth simulation study (Early Crossing with diminishing effect)

r15 <- c(.3, .2, .2, .1, .1, .15, .15, .2, .2, .2, .25)
r05 <- c(.2, .25, .25, .25, .25, .2, .2, .2, .2, .2, .2)
int5 <- c(rep(1/2, 10), Inf)

ff <- function(x) {
  ans <- Spwexp_point(x, rate=r15, intervals=int5, xmax=8) -  Spwexp_point(x, rate=r05, intervals=int5, xmax=8)
  return(ans)
}
cross5 <- uniroot(ff, interval=c(.1, 7.9))$root
cross5.y <- Spwexp_point(cross5, rate=r15, intervals=int5, xmax=8)

#####################################################
### Sixth simulation study
### have two crossings, but the first is a very
##    weak crossing while the second is a strong one


r16 <- c(.25, .3, .35, .35, .275, .15, .15, .1, .05, .05, .05)
r06 <- c(.275, .25, .25, .25, .25, .2, .2, .2, .2, .2, .2)
int6 <- c(rep(1/2, 10), Inf)

ff <- function(x) {
  ans <- Spwexp_point(x, rate=r16, intervals=int6, xmax=8) -  Spwexp_point(x, rate=r06, intervals=int6, xmax=8)
  return(ans)
}
cross61 <- uniroot(ff, interval=c(.1, 2))$root
cross62 <- uniroot(ff, interval=c(2.1, 7))$root
cross61.y <- Spwexp_point(cross61, rate=r16, intervals=int6, xmax=8)
cross62.y <- Spwexp_point(cross62, rate=r16, intervals=int6, xmax=8)


######################################################

x.max <- 8
s01 <- Spwexp(rate=r01, intervals=int1, xmax=x.max) 
s11 <- Spwexp(rate=r11, intervals=int1, xmax=x.max) 
s02 <- Spwexp(rate=r02, intervals=int2, xmax=x.max) 
s12 <- Spwexp(rate=r12, intervals=int2, xmax=x.max)
s03 <- Spwexp(rate=r03, intervals=int3, xmax=x.max) 
s13 <- Spwexp(rate=r13, intervals=int3, xmax=x.max) 
s04 <- Spwexp(rate=r04, intervals=int4, xmax=x.max) 
s14 <- Spwexp(rate=r14, intervals=int4, xmax=x.max) 
s05 <- Spwexp(rate=r05, intervals=int5, xmax=x.max) 
s15 <- Spwexp(rate=r15, intervals=int5, xmax=x.max) 
s06 <- Spwexp(rate=r06, intervals=int6, xmax=x.max) 
s16 <- Spwexp(rate=r16, intervals=int6, xmax=x.max) 

line.width <- 2
xleft.pos <- -1.2
#postscript("~/PWExpSimulationScenarios.eps", horizontal=FALSE, width=6.5, height=5.5)
par(mfrow=c(3,2), mar=c(2.6, 3.0, .25, .25))
plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s01$x, s01$y, col="black", lwd=line.width)
lines(s11$x, s11$y, col="red", lwd=line.width)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(a)   Scenario 1   ", bty='n', cex=1.3)

plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s02$x, s02$y, col="black", lwd=line.width)
lines(s12$x, s12$y, col="red", lwd=line.width)
arrows(cross2, cross2.y + .3, cross2, cross2.y, length=.15, lwd=2)
text(cross2 - .05, cross2.y + .37, "Cross", cex=1.2)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(b)   Scenario 2   ", bty='n', cex=1.3)


plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s03$x, s03$y, col="black", lwd=line.width)
lines(s13$x, s13$y, col="red", lwd=line.width)
arrows(cross3, cross3.y - .3, cross3, cross3.y, length=.15, lwd=2)
text(cross3 - .05, cross3.y - .35, "Cross", cex=1.2)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(c)   Scenario 3   ", bty='n', cex=1.3)


plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s04$x, s04$y, col="black", lwd=line.width)
lines(s14$x, s14$y, col="red", lwd=line.width)
arrows(cross4, cross4.y - .3, cross4, cross4.y, length=.15, lwd=2)
text(cross4 - .05, cross4.y - .35, "Cross", cex=1.2)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(d)   Scenario 4   ", bty='n', cex=1.3)


plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s05$x, s05$y, col="black", lwd=line.width)
lines(s15$x, s15$y, col="red", lwd=line.width)
arrows(cross5, cross5.y - .3, cross5, cross5.y, length=.15, lwd=2)
text(cross5 - .05, cross5.y - .35, "Cross", cex=1.2)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(e)   Scenario 5   ", bty='n', cex=1.3)


plot(0,0, type="n", xlim=c(0,x.max), ylim=c(0, 1), las=1, xlab="", ylab="",
     xaxt='n', cex.axis=0.9, yaxt='n')
lines(s06$x, s06$y, col="black", lwd=line.width)
lines(s16$x, s16$y, col="red", lwd=line.width)
arrows(cross61, cross61.y - .3, cross61, cross61.y, length=.15, lwd=2)
arrows(cross62, cross62.y + .3, cross62, cross62.y, length=.15, lwd=2)
text(cross61 - .05, cross61.y - .35, "Cross", cex=1.2)
text(cross62 - .05, cross62.y + .37, "Cross", cex=1.2)
text(x=4, y=-.23, labels="Time", xpd=NA, cex=1.1)
axis(1, at=c(0, 2, 4, 6, 8), labels=TRUE, tick=TRUE, padj=0, tcl=1, tck=.03, mgp=c(3,0.1,0), cex.axis=1)
text(x=xleft.pos, y=.5, labels="Survival Probability", xpd=NA, cex=1.1, srt=90)
axis(2, at=c(0, .2, .4, .6, .8, 1), labels=TRUE, padj=0, tcl=1, tck=-.03, mgp=c(3,0.5,0), cex.axis=1, las=1)
legend("topright", legend="(f)   Scenario 6   ", bty='n', cex=1.3)

#dev.off()






