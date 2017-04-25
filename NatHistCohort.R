# CF McQuaid
# 21/04/2017
# Cohort model of TB including regression and a slow stream
########################################
#         ^     ^     ^     ^          #
#      -> Q <-> S <-> C <-> Y -> M     # Q - Quiescent, S - Subclinical, C - Minimal (C+S-), Y - Disseminated (C+S+, symptomatic), M - Death (due to disease)
#     /   ^     ^                      #
#    E    |     |                      # E - Infectious, T - cumulative number disseminated
#     \   v     v                      #
#      -> K <-> Z                      # K- Quiescent (slow stream), Z- Subclinical (slow stream),
#         v     v                      #
########################################
setwd("C:/Users/cfmcquaid/Simulations")
library("reshape2"); library("deSolve"); library("ggplot2"); #source("Completed/multiplot.R")
# Parameter sets: Average (), Wax (regression), Slow (slow), Full (regression, slow)
# Ax - rate from compartment A to compartment X, Omega - background mortality
paramA <- c(Eq=1.00, Ek=0.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0.02)
paramW <- c(Eq=1.00, Ek=0.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0.02)
paramS <- c(Eq=0.10, Ek=0.90, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.05, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0.02)
paramF <- c(Eq=0.10, Ek=0.90, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.25, Zs=0.05, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0.02)
# Initial states of compartments
state <- c(E=100000, Q=0, K=0, S=0, Z=0, C=0, Y=0, T=0, M=0)
# Timespan for simulation
times <- seq(0, 25, by = 1)
# ODEfunction
regr <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
   # rates of change
   dE <- - (Eq + Ek)*E
   dQ <- + Eq*E + Sq*S + Kq*K - (Qs + Qk + Omega)*Q
   dK <- + Ek*E + Zk*Z + Qk*Q - (Kz + Kq + Omega)*K
   dS <- + Qs*Q + Cs*C + Zs*Z - (Sc + Sq + Sz + Omega)*S
   dZ <- + Kz*K + Sz*S - (Zk + Zs + Omega)*Z
   dC <- + Sc*S + Yc*Y - (Cs + Cy + Omega)*C
   dY <- + Cy*C -(Yc + Ym + Omega)*Y
   dT <- + Cy*C
   dM <- + Ym*Y
   # return the rate of change
   list(c(dE, dQ, dK, dS, dZ, dC, dY, dT, dM))
  })
}
# Calculations
outA <- ode(y = state, times = times, func = regr, parms = paramA); outW <- ode(y = state, times = times, func = regr, parms = paramW); outS <- ode(y = state, times = times, func = regr, parms = paramS); outF <- ode(y = state, times = times, func = regr, parms = paramF)
outA <- as.data.frame(outA); outW <- as.data.frame(outW); outS <- as.data.frame(outS); outF <- as.data.frame(outF)
# Formatting the data for plotting - putting into a melted data.fame, with additional columns for the source matrix and the incidence rate
outA$rate <- outA$S*paramA["Sc"]; outW$rate <- outW$S*paramW["Sc"]; outS$rate <- outS$S*paramS["Sc"]; outF$rate <- outF$S*paramF["Sc"];
outA <- melt(outA, id.vars = c("time")); outW <- melt(outW, id.vars = c("time")); outS <- melt(outS, id.vars = c("time")); outF <- melt(outF, id.vars = c("time"))
outA$source <- "A"; outW$source <- "W"; outS$source <- "S"; outF$source <- "F"
out <- rbind(outA, outW, outS, outF)
# Plot output
ggplot(out[out$variable %in%c("rate"), ], aes(time, value)) + geom_point(size=2) + labs(x = "Time", y = "Incidence") + facet_grid(source ~ . , scales = "free")
# Fitting parameters and state: output (shown here for the full parameter set) for IDing compartment sizes to generate reasonable initial state, then calculating values so that 2.5-5% of individuals are diseased or dead (due to TB) after 5 years
Qsz <- out[which(out$source == "F" & out$variable == "Q" & out$time == 5), ]; Ssz <- out[which(out$source == "F" & out$variable == "S" & out$time == 5), ]; Ksz <- out[which(out$source == "F" & out$variable == "K" & out$time == 5), ]; Zsz <- out[which(out$source == "F" & out$variable == "Z" & out$time == 5), ]; Csz <- out[which(out$source == "F" & out$variable == "C" & out$time == 5), ]; Ysz <- out[which(out$source == "F" & out$variable == "Y" & out$time == 5), ]; Msz <- out[which(out$source == "F" & out$variable == "M" & out$time == 5), ];
(Csz$value + Ysz$value + Msz$value) / (Qsz$value + Ssz$value + Ksz$value + Zsz$value + Csz$value + Ysz$value + Msz$value)
