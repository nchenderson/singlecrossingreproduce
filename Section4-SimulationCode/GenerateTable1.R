
## Code used to generate Table 1 from main manuscript
load("~/Section4-SimulationResults/run1_n200.RData")
load("~/Section4-SimulationResults/run2_n200.RData")
load("~/Section4-SimulationResults/run3_n200.RData")
load("~/Section4-SimulationResults/run4_n200.RData")
load("~/Section4-SimulationResults/run5_n200.RData")
load("~/Section4-SimulationResults/run6_n200.RData")

RR <- matrix(0, 18, 9)
k <- 0
RR[k + 1,] <- c(Results1[2,4], ResultsKM1[1,4], Results1[3,4], ResultsKM1[2,4], Results1[4,4], 
                ResultsKM1[3,4], Results1[1,4], Results1[5,4], Results1[6,4])
RR[k + 2,] <- c(Results2[2,4], ResultsKM2[1,4], Results2[3,4], ResultsKM2[2,4], Results2[4,4], 
                ResultsKM2[3,4], Results2[1,4], Results2[5,4], Results2[6,4])
RR[k + 3,] <- c(Results3[2,4], ResultsKM3[1,4], Results3[3,4], ResultsKM3[2,4], Results3[4,4], 
                ResultsKM3[3,4], Results3[1,4], Results3[5,4], Results3[6,4])
RR[k + 4,] <- c(Results4[2,4], ResultsKM4[1,4], Results4[3,4], ResultsKM4[2,4], Results4[4,4], 
                ResultsKM4[3,4], Results4[1,4], Results4[5,4], Results4[6,4])
RR[k + 5,] <- c(Results5[2,4], ResultsKM5[1,4], Results5[3,4], ResultsKM5[2,4], Results5[4,4], 
                ResultsKM5[3,4], Results5[1,4], Results5[5,4], Results5[6,4])
RR[k + 6,] <- c(Results6[2,4], ResultsKM6[1,4], Results6[3,4], ResultsKM6[2,4], Results6[4,4], 
                ResultsKM6[3,4], Results6[1,4], Results6[5,4], Results6[6,4])

rm(list=c("Results1", "ResultsKM1", "Results2", "ResultsKM2", "Results3", "ResultsKM3",
          "Results4", "ResultsKM4", "Results5", "ResultsKM5", "Results6", "ResultsKM6"))
load("~/Section4-SimulationResults/run1.RData")
load("~/Section4-SimulationResults/run2.RData")
load("~/Section4-SimulationResults/run3.RData")
load("~/Section4-SimulationResults/run4.RData")
load("~/Section4-SimulationResults/run5.RData")
load("~/Section4-SimulationResults/run6.RData")
k <- 6
RR[k + 1,] <- c(Results1[2,4], ResultsKM1[1,4], Results1[3,4], ResultsKM1[2,4], Results1[4,4], 
                ResultsKM1[3,4], Results1[1,4], Results1[5,4], Results1[6,4])
RR[k + 2,] <- c(Results2[2,4], ResultsKM2[1,4], Results2[3,4], ResultsKM2[2,4], Results2[4,4], 
                ResultsKM2[3,4], Results2[1,4], Results2[5,4], Results2[6,4])
RR[k + 3,] <- c(Results3[2,4], ResultsKM3[1,4], Results3[3,4], ResultsKM3[2,4], Results3[4,4], 
                ResultsKM3[3,4], Results3[1,4], Results3[5,4], Results3[6,4])
RR[k + 4,] <- c(Results4[2,4], ResultsKM4[1,4], Results4[3,4], ResultsKM4[2,4], Results4[4,4], 
                ResultsKM4[3,4], Results4[1,4], Results4[5,4], Results4[6,4])
RR[k + 5,] <- c(Results5[2,4], ResultsKM5[1,4], Results5[3,4], ResultsKM5[2,4], Results5[4,4], 
                ResultsKM5[3,4], Results5[1,4], Results5[5,4], Results5[6,4])
RR[k + 6,] <- c(Results6[2,4], ResultsKM6[1,4], Results6[3,4], ResultsKM6[2,4], Results6[4,4], 
                ResultsKM6[3,4], Results6[1,4], Results6[5,4], Results6[6,4])
rm(list=c("Results1", "ResultsKM1", "Results2", "ResultsKM2", "Results3", "ResultsKM3",
          "Results4", "ResultsKM4", "Results5", "ResultsKM5", "Results6", "ResultsKM6"))

load("~/Section4-SimulationResults/run1_n800.RData")
load("~/Section4-SimulationResults/run2_n800.RData")
load("~/Section4-SimulationResults/run3_n800.RData")
load("~/Section4-SimulationResults/run4_n800.RData")
load("~/Section4-SimulationResults/run5_n800.RData")
load("~/Section4-SimulationResults/run6_n800.RData")
k <- 12
RR[k + 1,] <- c(Results1[2,4], ResultsKM1[1,4], Results1[3,4], ResultsKM1[2,4], Results1[4,4], 
                ResultsKM1[3,4], Results1[1,4], Results1[5,4], Results1[6,4])
RR[k + 2,] <- c(Results2[2,4], ResultsKM2[1,4], Results2[3,4], ResultsKM2[2,4], Results2[4,4], 
                ResultsKM2[3,4], Results2[1,4], Results2[5,4], Results2[6,4])
RR[k + 3,] <- c(Results3[2,4], ResultsKM3[1,4], Results3[3,4], ResultsKM3[2,4], Results3[4,4], 
                ResultsKM3[3,4], Results3[1,4], Results3[5,4], Results3[6,4])
RR[k + 4,] <- c(Results4[2,4], ResultsKM4[1,4], Results4[3,4], ResultsKM4[2,4], Results4[4,4], 
                ResultsKM4[3,4], Results4[1,4], Results4[5,4], Results4[6,4])
RR[k + 5,] <- c(Results5[2,4], ResultsKM5[1,4], Results5[3,4], ResultsKM5[2,4], Results5[4,4], 
                ResultsKM5[3,4], Results5[1,4], Results5[5,4], Results5[6,4])
RR[k + 6,] <- c(Results6[2,4], ResultsKM6[1,4], Results6[3,4], ResultsKM6[2,4], Results6[4,4], 
                ResultsKM6[3,4], Results6[1,4], Results6[5,4], Results6[6,4])

round(RR, 4)

library(xtable)
Cnames <- c("RMST Diff", "RMST Diff (KM)", "S2 Diff", "S2 Diff(KM)", "S4 Diff", "S4 Diff(KM)",
            "Cross Time", "RRML Diff", "S(theta)")
cap <- "Mean-squared error for different quantities using the single-crossing constrained estimator
        and the Kaplan-Meier estimator. "

colnames(RR) <- Cnames
xtable(RR[,c(1:7,9,8)], digits=rep(4, 10), caption=cap)




