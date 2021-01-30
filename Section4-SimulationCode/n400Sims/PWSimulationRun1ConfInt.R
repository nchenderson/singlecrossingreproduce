library(boot)
library(DelayedSurvFit)

rpwexp <- function(n, rate, intervals){
  
  n.int <- length(rate)
  if(n.int==1){
    return(rexp(n, rate))
  }
  tx <- 0
  j <- 1
  times <- array(0,n)
  timex <- cumsum(intervals)
  indx <- rep(TRUE,n)
  for(i in 1:n.int){
    nindx <- sum(indx)
    if (nindx==0) {
      break
    } 
    increment <- rexp(nindx,rate[i])
    times[indx] <- tx + increment
    if (i < n.int){
      tx <- timex[i]
      indx <- (times > timex[i])
    }
  }
  return(times)
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

theta.fun <- function(obj, trt) {
  ## Need to output the following quantities
  # (1) theta
  # (2) RMST
  # (3) Survival Prob difference
  # (4) S_{a}(\theta)
  # (5) RRML(\theta)
  obj.ds <- DelayedSurvFit(obj[,1], obj[,2], trt)
  ans <- rep(0, 6)
  ans[1] <- ifelse(obj.ds$theta <= min(obj[,1]), 0, obj.ds$theta)
  #ans[2] <- rnorm(1)
  rmst.estim <- rmst(obj.ds, tau0=5, tau1=5)
  rrml.estim <- rrml(obj.ds, t0=ds.obj$theta, t1=ds.obj$theta, tau0=5, tau1=5)
  SS0 <- SurvFn(obj.ds, arm=0)
  SS1 <- SurvFn(obj.ds, arm=1)
  ans[2] <- rmst.estim$rmst1 - rmst.estim$rmst0
  ans[3] <- SS1(1) - SS0(1)
  ans[4] <- SS1(3) - SS0(3)
  ans[5] <- median(SS1(ans[1]), SS0(ans[1]))
  ans[6] <- rrml.estim$rrml0
  return(ans)
}

############################################################
# Simulation Study 1 (one survival curve dominates the other)

library(DelayedSurvFit)
library(survival)
#r01 <- c(.8, .5, .2, .1, .05)
r01 <- c(.7, .5, .2, .1, .05)
r11 <- c(.5, .3, .15, .075, .025)
int1 <- c(rep(1, 4), Inf)

rmst0.true <- integrate(Spwexp, lower=0, upper=5, rate=r01, intervals=int1, subdivisions=1000L)$value
rmst1.true <- integrate(Spwexp, lower=0, upper=5, rate=r11, intervals=int1, subdivisions=1000L)$value
survfdiff.true <- Spwexp(1, rate=r11, intervals=int1) - Spwexp(1, rate=r01, intervals=int1) 
survsdiff.true <- Spwexp(3, rate=r11, intervals=int1) - Spwexp(3, rate=r01, intervals=int1)
survcross.true <- 1
rmstdiff.true <- rmst1.true - rmst0.true
## RRML.true <- 

n <- 50
nreps <- 1
n.bootstrap <- 2  ## 100 works pretty well

set.seed(57216)
theta.hat <- min.obs <- gamma.hat <- RMST0 <- RMST1 <- RMST.KM0 <- RMST.KM1 <- RMSTDiff <- S2Diff <- S4Diff <- Stheta <- rep(0, nreps)
CI.theta <- CI.RMSTDiff <- CI.SF <- CI.SS <- CI.Stheta <- matrix(0, nrow=nreps, ncol=2)
for(k in 1:nreps) {
  times0.true <- rpwexp(n, rate=r01, intervals=int1)
  times1.true <- rpwexp(n, rate=r11, intervals=int1)
  cens0 <- runif(n, min=4, max=6)
  cens1 <- runif(n, min=4, max=6)
  event0 <- as.numeric(times0.true < cens0)
  event1 <- as.numeric(times1.true < cens1)
  times0 <- pmin(times0.true, cens0)
  times1 <- pmin(times1.true, cens1)
  trt <- rep(c(0, 1), each=n)
  
  times <- c(times0, times1)
  events <- c(event0, event1)
  
  ds.obj <- DelayedSurvFit(times, events, trt)
  ds.obj$theta <- ifelse(ds.obj$theta <= min(times), 0, ds.obj$theta)
  
  ## Record results
  theta.hat[k] <- ds.obj$theta
  min.obs[k] <- min(times)
  gamma.hat[k] <- ds.obj$gamma
  tmp.rmst <- rmst(ds.obj, tau0=5, tau1=5)
  RMST0[k] <- tmp.rmst$rmst0
  RMST1[k] <- tmp.rmst$rmst1
  RMSTDiff[k] <- tmp.rmst$rmst1 - tmp.rmst$rmst0
  SS0 <- SurvFn(ds.obj, arm=0)
  SS1 <- SurvFn(ds.obj, arm=1)
  S2Diff[k] <- SS1(1) - SS0(1)
  S4Diff[k] <- SS1(3) - SS0(3)
  Stheta[k] <- median(SS1(ds.obj$theta), SS0(ds.obj$theta))
  
  ## Boostrap and confidence intervals
  #dat <- cbind(times, events)
  #bt.info <- censboot(dat, theta.fun, R = n.bootstrap, trt = trt)
  ## bt.info will be a matrix with dims nreps x n.statistics
  
  #CI.theta[k,] <- quantile(as.vector(bt.info$t[,1]), probs=c(0.025, 0.975))
  #CI.RMSTDiff[k,] <- quantile(as.vector(bt.info$t[,2]), probs=c(0.025, 0.975))
  #CI.SF[k,] <- quantile(as.vector(bt.info$t[,3]), probs=c(0.025, 0.975))
  #CI.SS[k,] <- quantile(as.vector(bt.info$t[,4]), probs=c(0.025, 0.975))
  #CI.Stheta[k,] <- quantile(as.vector(bt.info$t[,5]), probs=c(0.025, 0.975))
  
  
  ## Compute KM estimates
  # aa0 <- survfit(Surv(times[1:n],events[1:n])~1)
  #  aa1 <- survfit(Surv(times[(n+1):(2*n)], events[(n+1):(2*n)])~1)
  #  tmpp.rmst <- rmst2(times, events, arm=trt, tau=5)
  #  RMST.KM0[k] <- tmpp.rmst$RMST.arm0$rmst[1]
  #  RMST.KM1[k] <- tmpp.rmst$RMST.arm1$rmst[1]
  # tst <- RMST.KM0[k]/RMST0[k] > 1.1 | RMST.KM0[k]/RMST0[k] < 0.9
  #  tst <- ds.obj$theta > 1
  #  tst <- TRUE
  #  if(tst) {
  
  #      plot(aa0, conf.int=FALSE, main=k)
  #      lines(aa1$time, aa1$surv, col="red", type="s")
  #      lines(ds.obj$times0, ds.obj$surv0, lty=2, type="s")
  #      lines(ds.obj$times1, ds.obj$surv1, lty=2, col="red", type="s")
  #    
  #  }
  print(k)
  
}
## Want most estimates of theta to be less than 0.1 or 0.2





