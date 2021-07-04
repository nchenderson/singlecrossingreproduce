library(DelayedSurvFit)
data("NivoIpili")
head(NivoIpili)
ds.obj <- DelayedSurvFit(times=NivoIpili$times, events=NivoIpili$evnt, 
                         trt=NivoIpili$trt, inner.iter=100)


S1.hat <- SurvFn(ds.obj, arm=1)
S0.hat <- SurvFn(ds.obj, arm=0)

s.theta.orig <- (S1.hat(ds.obj$theta) + S0.hat(ds.obj$theta))/2
tmp <- rmst(ds.obj, tau0=36, tau1=36)
rmst.diff.orig <- tmp$rmst1 - tmp$rmst0

## s.theta.orig = 0.72876
## rmst.diff.orig = 1.4826552

## Now do bootstrap
R <- 200
p.star <- 0.6 ## Think about using p.star = 0.6 or p.star = 0.5 instead.
eta <- rmd <- sthet <- rep(NA, R)
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
                           trt=TmpDat$trt, max.times=100, inner.iter=50)
  
  SB1.hat <- SurvFn(bt.obj, arm=1)
  SB0.hat <- SurvFn(bt.obj, arm=0)
  s.theta <- (SB1.hat(bt.obj$theta) + SB0.hat(bt.obj$theta))/2
  tmp <- rmst(bt.obj, tau0=36, tau1=36)
  rmst.diff <- tmp$rmst1 - tmp$rmst0
  eta[k] <- min(rmst.diff, s.theta - p.star)
  rmd[k] <- rmst.diff 
  sthet[k] <- s.theta
  cat("Iteration", k, "\n")
}


quantile(eta, probs=0.05)

eta2 <- pmin(rmd, sthet - 0.6)
round(quantile(eta2, probs=0.05), 3)

summary(rmd)
summary(sthet - 0.6)
round(min(rmst.diff.orig, s.theta.orig - p.star), 3)



