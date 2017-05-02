# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version 26/04/2017
# Cohort model of TB including regression and a slow stream
########################################
#         ^     ^     ^     ^          #
#      -> Q <-> S <-> C <-> Y -> M     # Q - Quiescent, S - Subclinical, C - Minimal (C+S-), Y - Disseminated (C+S+, symptomatic), M - Death (due to disease)
#     /   ^     ^                      #
#    E    |     |                      # E - Recently Infected, T - Cumulative incidence
#     \   v     v                      #
#      -> K <-> Z                      # K - Quiescent (slow stream), Z- Subclinical (slow stream),
#         v     v                      #
########################################
##set wd depending on who is coding/playing
setwd("C:/Users/cfmcquaid/Simulations")
#setwd("/Users/ReinHouben/Filr/ifolder/Applications (funding and jobs)/2016 - ERC starting grant/Rmodel_nathis/NatHistCoh")

#Next steps
# 1. get the prop inc disease after 5 years to 5%
###   How: play with Eq and Ek values, OR increase progression parameters
# 2. achieve 1 + start higher and sharper decrease in early years
###   How: change initial states to more realistic point, based on what prev of disease would be at screening
###   Also increase progression rates, while decreasing proportion that goes into the fast pathway

library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr");
# Parameter sets: Average (), Wax (regression), Slow (slow), Full (regression, slow)
# Ax - rate from compartment A to compartment X, Omega - background mortality
paramA <- c(Eq=1.00, Ek=0.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0)
paramW <- c(Eq=1.00, Ek=0.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0)
paramS <- c(Eq=0.10, Ek=0.90, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.05, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0)
paramF <- c(Eq=0.10, Ek=0.90, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.25, Zs=0.05, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0)
# Initial states of compartments
state <- c(E=100000, Q=0, K=0, S=0, Z=0, C=0, Y=0, M=0, T=0)
# Timespan for simulation
times <- seq(0, 25, by = 1)
# Timespan for burn
timeb <- 5
# ODE function
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
   dM <- + Ym*Y
   dT <- + Sc*S / (E + Q + K + S + Z + C + Y + M)
   # return the rate of change
   list(c(dE, dQ, dK, dS, dZ, dC, dY, dM, dT))
  })
}
# Burn off the first X years to establish the proportion of individuals in each compartment, removing diseased
burn <- function(tb, state, fxn, parameters){
  out <- ode(y = state, times = seq(0, tb, by = 1), func = fxn, parms = parameters)
  stateb <- c(out[tb, 2], out[tb, 3], out[tb, 4], out[tb, 5], out[tb, 6], C=0, Y=0, M=0, T=0)
  return(stateb)
}
# Calculation function
calc <- function(t, tb, state, fxn, parameters, source){
  stateb <- burn(tb = tb, state = state, fxn = regr, parameters = parameters)
  out <- ode(y = stateb, times = t, func = fxn, parms = parameters)
  out <- as.data.frame(out)
  # Calculating incidence by subtracting the previous year's cumulative incidence from this year's
  out$inc <- head(c(out$T,0) - c(0,out$T),-1)
  # Removing the first data point at time zero as we only sample after a year
  out <- out[-c(1), ]
  # Formatting the data for plotting - putting into a melted data.fame, with additional columns for the source matrix
  out <- melt(out, id.vars = c("time"))
  out$source <- source
  return(out)
}
outA <- calc(t = times, tb = timeb, state = state, fxn = regr, parameters = paramA, source = "A")
outW <- calc(t = times, tb = timeb, state = state, fxn = regr, parameters = paramW, source = "W")
outS <- calc(t = times, tb = timeb, state = state, fxn = regr, parameters = paramS, source = "S")
outF <- calc(t = times, tb = timeb, state = state, fxn = regr, parameters = paramF, source = "F")
out <- rbind(outA, outW, outS, outF)
# Plot output
theme_set(theme_bw())
##all scenarios
ggplot(out[out$variable %in% c("inc"), ], aes(time, value)) + geom_point(size=2) + labs(x = "Time", y = "Incidence") + facet_grid(source ~ . , scales = "free")
##slow and slow+regression only (S+F)
out2 <- subset(out, source=='F' | source == 'S', select=time:source)  
#ggplot(out2[out2$variable %in% c("inc"), ], aes(time, value)) + geom_point(size=2) + labs(x = "Time", y = "Incidence") + facet_grid(source ~ . , scales = "free")

# Fitting parameters and state: proportion of individuals disease or dead (due to TB) after Otime years, using parameter set Osource
Osource <- "F"; Otime <- 5
Osize <- out[which(out$source == Osource & out$variable %in% c("E", "Q", "S", "K", "Z",  "C", "Y", "M") & out$time == Otime), ]
Oinc <- out[which(out$source == Osource & out$variable %in% c("T") & out$time == Otime), ]
# Point prevalence: total incidence excluding those that may have regressed
Pp <- sum(Osize$value[Osize$variable %in% c("C", "Y", "M")]) / sum(Osize$value)
# Cumulative incidence: total incidence including those that may have regressed
Ci <- Oinc$value

Pp
Ci

# Tornado plot
torn <- function(t, tb, tc, state, fxn, parameters){
  # Store data
  data <- rbind(parameters, parameters)
  rownames(data) <- c('+10%', '-10%')   
  # Compare to default data set
  stateb <- burn(tb = tb, state = state, fxn = regr, parameters = parameters)
  out <- ode(y = stateb, times = t, func = fxn, parms = parameters);
  def <- out [tc + 1, "T"]
    for (i in 1:16){
    # Increasing and decreasing each parameter in turn
    parametersM = parameters; parametersL = parameters
    parametersM[i] = parametersM[i] + 0.1*parametersM[i]; parametersL[i] = parametersL[i] - 0.1*parametersL[i]
    stateM <- burn(tb = tb, state = state, fxn = regr, parameters = parametersM);stateL <- burn(tb = tb, state = state, fxn = regr, parameters = parametersL)
    outM <- ode(y = stateM, times = t, func = fxn, parms = parametersM); outL <- ode(y = stateL, times = t, func = fxn, parms = parametersL)
    outM <- def - outM [tc + 1, "T"]; outL <- def - outL [tc + 1, "T"]
    data[1, i] <- outM; data[2, i] <- outL; 
  }
  return(data)
}
data <- torn(t = times, tb = timeb, tc = Otime, state = state, fxn = regr, parameters = paramF)
# For plotting '%' on x-axis
x <- seq(-0.01,0.01, length=10)
ORD = order(abs(data[2,] - data[1,]))
barplot(data[1,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('black'))
barplot(data[2,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('white'), add = TRUE)
axis(1, at=pretty(x), lab=paste0(pretty(x) * 100," %"), las=TRUE)
