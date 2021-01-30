library(DelayedSurvFit)
library(survival)


rpwexp <- function(n, rate, intervals)  {
  
  ff <- function(x, u) {
    Spwexp_point(x, rate=rate, intervals=intervals, xmax=500) - u
  }
  uu <- runif(n)
  yy <- rep(0, n)
  for(k in 1:n) {
    yy[k] <- uniroot(ff, lower=0, upper=500, u=uu[k])$root
  }
  return(yy)
}


Spwexp <- function(x, rate=1, intervals) {
  tau <- length(x)
  HH <- x
  for(k in 1:tau) {
    interval.idx <- max(which(x[k] > cumsum(c(0,intervals))))
    idx <- interval.idx - 1
    if(idx > 0) {
      HH[k] <- sum(rate[1:idx]*intervals[1:idx]) + (x[k] - sum(intervals[1:idx]))*rate[interval.idx]
    } else {
      HH[k] <- x[k]*rate[1]
    }
  }
  return(exp(-HH))
}


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

###################################################################
## Quantities to look at:
# (1) \theta
# (2) RMST
# (3) diff: S_{1}(2) - S_{0}(2)
# (4) diff: S_{1}(4) - S_{0}(4)
# (5) S_{a}( \theta)
# (6) RRML(\theta)
# (7) \gamma

## Measures to record
# (1) MSE
# (2) Coverage
# (3) proportion of sign errors

############################################################
# Simulation Study 1 (one survival curve dominates the other)

library(DelayedSurvFit)
library(survival)

r12 <- c(.3, .25, .2, .1, .05, .05)
r02 <- c(.1, .15, .25, .2, .2, .2)
int2 <- c(rep(1, 5), Inf)

ff <- function(x) {
  ans <- Spwexp(x, rate=r12, intervals=int2) -  Spwexp(x, rate=r02, intervals=int2)
  return(ans)
}
cross.point <- uniroot(ff, interval=c(.1, 7.9))$root

tau <- 7  
t1 <- 2
t2 <- 4
rmst0.true <- integrate(Spwexp, lower=0, upper=tau, rate=r02, intervals=int2, subdivisions=1000L)$value
rmst1.true <- integrate(Spwexp, lower=0, upper=tau, rate=r12, intervals=int2, subdivisions=1000L)$value
sfdiff.true <- Spwexp(t1, rate=r12, intervals=int2) - Spwexp(t1, rate=r02, intervals=int2) 
ssdiff.true <- Spwexp(t2, rate=r12, intervals=int2) - Spwexp(t2, rate=r02, intervals=int2)
survcross.true <- Spwexp(cross.point, rate=r12, intervals=int2)
rmstdiff.true <- rmst1.true - rmst0.true
num1 <- integrate(Spwexp, lower=cross.point, upper=tau, rate=r12, intervals=int2, subdivisions=1000L)$value  
num0 <- integrate(Spwexp, lower=cross.point, upper=tau, rate=r02, intervals=int2, subdivisions=1000L)$value  

rrmldiff.true <- num1/survcross.true - num0/survcross.true ## need to change this

n <- 100
nreps <- 200

set.seed(57216)
theta.hat <- min.obs <- gamma.hat <- RMST0 <- RMST1 <- RMSTDiff <- WKMStat <- rep(0, nreps)
SFDiff <- SSDiff <- Stheta <- RRML0 <- RRML1 <- RRMLDiff <- rep(0, nreps)
RMST.KM0 <- RMST.KM1 <- RMSTDiff.KM <- SFDiff.KM <- SSDiff.KM <- rep(0, nreps)
for(k in 1:nreps) {
  times0.true <- rpwexp(n, rate=r02, intervals=int2)
  times1.true <- rpwexp(n, rate=r12, intervals=int2)
  
  cens0 <- runif(n, min=4, max=8)
  cens1 <- runif(n, min=4, max=8)
  event0 <- as.numeric(times0.true < cens0)
  event1 <- as.numeric(times1.true < cens1)
  times0 <- pmin(times0.true, cens0)
  times1 <- pmin(times1.true, cens1)
  trt <- rep(c(0, 1), each=n)
  
  times <- c(times0, times1)
  events <- c(event0, event1)
  
  ds.obj <- DelayedSurvFit(times, events, trt, inner.iter=100)
  
  ## Record results
  theta.hat[k] <- ds.obj$theta
  min.obs[k] <- min(times)
  gamma.hat[k] <- ds.obj$gamma
  tmp.rmst <- rmst(ds.obj, tau0=tau, tau1=tau)
  RMST0[k] <- tmp.rmst$rmst0
  RMST1[k] <- tmp.rmst$rmst1
  RMSTDiff[k] <- tmp.rmst$rmst1 - tmp.rmst$rmst0
  SS0 <- SurvFn(ds.obj, arm=0)
  SS1 <- SurvFn(ds.obj, arm=1)
  SFDiff[k] <- SS1(t1) - SS0(t1)
  SSDiff[k] <- SS1(t2) - SS0(t2)
  Stheta[k] <- median(SS1(ds.obj$theta), SS0(ds.obj$theta))
  tmp.rrml <- rrml(ds.obj, t0=ds.obj$theta, t1=ds.obj$theta, tau0=tau, tau1=tau)
  RRML0[k] <- tmp.rrml$rrml0
  RRML1[k] <- tmp.rrml$rrml1
  RRMLDiff[k] <- tmp.rrml$rrml1 - tmp.rrml$rrml0
  
  ## Now compute the Kaplan-Meier estimate of each quantity except the crossing
  kmfit0 <- survfit(Surv(ds.obj$discretized.times[trt==0], events[trt==0]) ~ 1)
  kmfit1 <- survfit(Surv(ds.obj$discretized.times[trt==1], events[trt==1]) ~ 1)
  
  KM0 <- stepfun(kmfit0$time, c(1, kmfit0$surv), right=min(kmfit0$surv))
  KM1 <- stepfun(kmfit1$time, c(1, kmfit1$surv), right=min(kmfit1$surv))
  SFDiff.KM[k] <- KM1(t1) - KM0(t1)
  SSDiff.KM[k] <- KM1(t2) - KM0(t2)
  
  #tmp <- rmst2(ds.obj$discretized.times, events, arm=trt, tau=5)
  #tmp <- rmst2(times, events, arm=trt, tau=5)
  S0 <- c(1, kmfit0$surv[kmfit0$time <= tau])
  S1 <- c(1, kmfit1$surv[kmfit1$time <= tau])
  
  xx0 <- c(0, kmfit0$time[kmfit0$time <= tau], tau)
  xx1 <- c(0, kmfit1$time[kmfit1$time <= tau], tau)
  rmst0 <- sum(S0*diff(xx0))
  rmst1 <- sum(S1*diff(xx1))
  
  RMSTDiff.KM[k] <- rmst1 - rmst0 #tmp$RMST.arm1$rmst[1] - tmp$RMST.arm0$rmst[1]
  print(k)
}

nquants <- 6
nquants.KM <- 3
Results <- matrix(0, nrow=nquants, ncol=5)
colnames(Results) <- c("mean", "median", "sd", "MSE", "True Value")
rownames(Results) <- c("Cross Time", "RMST.Diff", "S1.Diff", "S2.Diff", "RRML.Diff","S(cross)")
Results[1,] <- c(mean(theta.hat), median(theta.hat), sd(theta.hat), mean((theta.hat - cross.point)^2), cross.point)
Results[2,] <- c(mean(RMSTDiff), median(RMSTDiff), sd(RMSTDiff), mean((RMSTDiff - rmstdiff.true)^2), rmstdiff.true)
Results[3,] <- c(mean(SFDiff), median(SFDiff), sd(SFDiff), mean((SFDiff - sfdiff.true)^2), sfdiff.true)
Results[4,] <- c(mean(SSDiff), median(SSDiff), sd(SSDiff), mean((SSDiff - ssdiff.true)^2), ssdiff.true)
Results[5,] <- c(mean(RRMLDiff), median(RRMLDiff), sd(RRMLDiff), mean((RRMLDiff - rrmldiff.true)^2), rrmldiff.true)
Results[6,] <- c(mean(Stheta), median(Stheta), sd(Stheta), mean((Stheta - survcross.true)^2), survcross.true)

round(Results, 4)

ResultsKM <- matrix(0, nrow=nquants.KM, ncol=5)
colnames(ResultsKM) <- c("mean", "median", "sd", "MSE", "True Value")
rownames(ResultsKM) <- c("RMST.Diff", "S1.Diff", "S2.Diff")
ResultsKM[1,] <- c(mean(RMSTDiff.KM), median(RMSTDiff.KM), sd(RMSTDiff.KM), mean((RMSTDiff.KM - rmstdiff.true)^2), rmstdiff.true)
ResultsKM[2,] <- c(mean(SFDiff.KM), median(SFDiff.KM), sd(SFDiff.KM), mean((SFDiff.KM - sfdiff.true)^2), sfdiff.true)
ResultsKM[3,] <- c(mean(SSDiff.KM), median(SSDiff.KM), sd(SSDiff.KM), mean((SSDiff.KM - ssdiff.true)^2), ssdiff.true)

round(ResultsKM, 4)

Results2 <- Results
ResultsKM2 <- ResultsKM

save(theta.hat, RMSTDiff, SFDiff, SSDiff, RRMLDiff, Stheta, RMSTDiff.KM, SFDiff.KM, SSDiff.KM,
     rmstdiff.true, sfdiff.true, ssdiff.true, rrmldiff.true, Results2, ResultsKM2, 
     file="~/Section4-SimulationResults/run2_n200.RData")







