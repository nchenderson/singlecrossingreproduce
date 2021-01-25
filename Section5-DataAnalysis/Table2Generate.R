#####################################################
####  Code to generate the numbers shown in Table 2
####  of the main manuscript
#####################################################

library(DelayedSurvFit)
library(xtable)
data("NivoIpiliReconstruct")
head(NivoIpili)
ds.obj <- DelayedSurvFit(times=NivoIpili$times, events=NivoIpili$evnt, 
                         trt=NivoIpili$trt, inner.iter=100)


S1.hat <- SurvFn(ds.obj, arm=1)
S0.hat <- SurvFn(ds.obj, arm=0)
Results <- matrix(0, nrow=10, ncol=1)
rownames(Results) <- c("theta", "S(theta)", "RMST(36) diff", "RRML(theta, 36) diff", 
                       "Delta S(6)", "Delta S(12)", "Delta S(24)",
                       "Delta S_c(12)", "Delta S_c(24)", "Delta S_c(36)")
Results[1,1] <- ds.obj$theta
s.theta <- (S1.hat(ds.obj$theta) + S0.hat(ds.obj$theta))/2
Results[2,1] <- s.theta
tmp <- rmst(ds.obj, tau0=36, tau1=36)
Results[3,1] <- tmp$rmst1 - tmp$rmst0
tmp <- rrml(ds.obj, t0=ds.obj$theta, t1=ds.obj$theta, tau0=36, tau1=36)
Results[4,1] <- tmp$rrml1 - tmp$rrml0 
Results[5,1] <- S1.hat(6) - S0.hat(6)
Results[6,1] <- S1.hat(12) - S0.hat(12)
Results[7,1] <- S1.hat(24) - S0.hat(24)
Results[8,1] <- S1.hat(12)/s.theta - S0.hat(12)/s.theta
Results[9,1] <- S1.hat(24)/s.theta - S0.hat(24)/s.theta
Results[10,1] <- S1.hat(36)/s.theta - S0.hat(36)/s.theta

## Now do bootstrap
R <- 10
BtValues <- matrix(0, nrow=10, ncol=R)
ind1 <- NivoIpili$trt == 1
ind0 <- NivoIpili$trt == 0
n1 <- sum(ind1)
n0 <- sum(ind0)
group1 <- which(ind1)
group0 <- which(ind0)

for(k in 1:R) {
  ss1 <- sample(group1, size=n1, replace=TRUE)
  ss0 <- sample(group0, size=n0, replace=TRUE)
  times1 <- NivoIpili$times[ss1]
  evt1 <- NivoIpili$evnt[ss1]
  times0 <- NivoIpili$times[ss0]
  evt0 <- NivoIpili$evnt[ss0]
  
  TmpDat <- data.frame(times=c(times1, times0), evnt=c(evt1, evt0), trt=c(rep(1, n1), rep(0, n0)))
  bt.obj <- DelayedSurvFit(times=TmpDat$times, events=TmpDat$evnt, 
                           trt=TmpDat$trt, inner.iter=50)
  
  SB1.hat <- SurvFn(bt.obj, arm=1)
  SB0.hat <- SurvFn(bt.obj, arm=0)
  BtValues[1,k] <- bt.obj$theta
  s.theta <- (SB1.hat(bt.obj$theta) + SB0.hat(bt.obj$theta))/2
  BtValues[2,k] <- s.theta
  tmp <- rmst(bt.obj, tau0=36, tau1=36)
  BtValues[3,k] <- tmp$rmst1 - tmp$rmst0
  tmp <- rrml(bt.obj, t0=bt.obj$theta, t1=bt.obj$theta, tau0=36, tau1=36)
  BtValues[4,k] <- tmp$rrml1 - tmp$rrml0 
  BtValues[5,k] <- SB1.hat(6) - SB0.hat(6)
  BtValues[6,k] <- SB1.hat(12) - SB0.hat(12)
  BtValues[7,k] <- SB1.hat(24) - SB0.hat(24)
  BtValues[8,k] <- SB1.hat(12)/s.theta - SB0.hat(12)/s.theta
  BtValues[9,k] <- SB1.hat(24)/s.theta - SB0.hat(24)/s.theta
  BtValues[10,k] <- SB1.hat(36)/s.theta - SB0.hat(36)/s.theta
  cat("Iteration", k, "\n")
}

UncertaintyResults <- t(apply(BtValues, 1, function(x) quantile(x, probs=c(0.025, 0.975))))

## Generate Latex table
library(xtable)
A <- cbind(Results, UncertaintyResults)
xtable(A, digits=rep(3, 4))


