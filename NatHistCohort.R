# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version 26/04/2017
# Cohort model of TB including regression and a slow stream
########################################
#         ^     ^     ^     ^          #
#      -> Q <-> S <-> C <-> Y -> M     # Q - Quiescent, S - Subclinical, C - Minimal (C+S-), Y - Disseminated (C+S+, symptomatic), M - Death (due to disease)
#     /   ^     ^                      #
#    E    |     |                      # E - Recently Infected, T - Cumulative incidence (using initial disease-free population)
#     \   v     v                      #
#      -> K <-> Z                      # K - Quiescent (slow stream), Z- Subclinical (slow stream),
#         v     v                      #
########################################
##set wd depending on who is coding/playing
#setwd("C:/Users/cfmcquaid/Simulations")
setwd("/Users/ReinHouben/Filr/ifolder/Applications (funding and jobs)/2016 - ERC starting grant/Rmodel_nathis/NatHistCoh")

### URL for how to upload changes: http://r-bio.github.io/intro-git-rstudio/
### to pull in Finn's version - type in shell: git pull upstream master
### to push my version - create pull request to Finn:  
#       commit (remember to tick box), then push
#       goto Github, pull requests, create pull request, if merge can be done on auto, just complete
#       Alternatively Finn complete pull request (he does not get email notification)

#Next steps
# 1. get the prop inc disease after 5 years to 5%
###   How: play with Eq values, OR increase progression parameters
# 2. achieve 1 + start higher and sharper decrease in early years
###   How: change initial states to more realistic point, based on what prev of disease would be at screening
###   Also increase progression rates, while decreasing proportion that goes into the fast pathway
# OTHER NEXT STEPS - 8th May 2017
##  ?can we calculate in model how often (on average) people would come in and out of states over course of say 2 years?
##  ?Can we calculate average time in compartment in each step visit? Should be simple competing risk model right?

library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr");
# Parameter sets: Average (), Wax (regression), Slow (slow), Full (regression, slow)
# Ax - rate from compartment A to compartment X, Omega - background mortality
paramA <- c(Eq=1.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0)
paramW <- c(Eq=1.00, Qs=0.03, Qk=0.00, Kz=0.00, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.00, Zs=0.00, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0)
paramS <- c(Eq=0.10, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.22, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.05, Cy=0.22, Cs=0.00, Ym=0.25, Yc=0.00, Omega=0)
paramF <- c(Eq=0.10, Qs=0.50, Qk=0.00, Kz=0.01, Kq=0.00, Sc=0.50, Sq=0.25, Sz=0.00, Zk=0.25, Zs=0.05, Cy=0.50, Cs=0.25, Ym=0.50, Yc=0.25, Omega=0)
# Initial states of compartments
state <- c(E=100000, Q=0, K=0, S=0, Z=0, C=0, Y=0, M=0, T=0)
# Timespan for simulation
times <- seq(0, 10, by = 1)
# Timespan for burn
timeb <- 2
# Timespan for tornado
timet <- 5
# ODE function
regr <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
   # rates of change
   dE <- - E
   dQ <- + Eq*E + Sq*S + Kq*K - (Qs + Qk + Omega)*Q
   dK <- + (1 - Eq)*E + Zk*Z + Qk*Q - (Kz + Kq + Omega)*K
   dS <- + Qs*Q + Cs*C + Zs*Z - (Sc + Sq + Sz + Omega)*S
   dZ <- + Kz*K + Sz*S - (Zk + Zs + Omega)*Z
   dC <- + Sc*S + Yc*Y - (Cs + Cy + Omega)*C
   dY <- + Cy*C -(Yc + Ym + Omega)*Y
   dM <- + Ym*Y
   dT <- + Sc*S / sum(state[1:5])
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
calc <- function(ts, tb, state, fxn, parameters, source){
  stateb <- burn(tb = tb, state = state, fxn = regr, parameters = parameters)
  out <- ode(y = stateb, times = ts, func = fxn, parms = parameters)
  out <- as.data.frame(out)
  # Calculating incidence
  out$inc <- parameters["Sc"]*out$S / (out$E + out$Q + out$S + out$K + out$Z)
  out$inc <- head(c(0,out$inc),-1) # We later remove the first element as it is zero if calculated using C&Y, but will be nonnegative if looking at the rate. Shift the index accordingly
  # Calculating incidence rate (https://en.wikipedia.org/wiki/Incidence_(epidemiology)), proportion of population that result in cases per year = (total # cases incl mort) / {[(C+Y)_t1-(C+Y)_t0]*(t1+t0)/2 + t1*(E+Q+S+K+Z)_t1}
  out$rate <- (out$C + out$Y + out$M) / (head(c(out$C,0) + c(out$Y,0) - c(0,out$C) - c(0,out$Y), -1)*head(c(out$time,0) - c(0,out$time), -1) / 2 + (out$E + out$Q + out$S + out$K + out$Z)*out$time)
  # Calculating the "relative risk" (see Vynnycky & Fine, 1997)
  # out$risk <- head(c(out$C,0) + c(out$Y,0) - c(0,out$C) - c(0,out$Y), -1)/(out$C[2] + out$Y[2])
  # out$risk <- parameters["Sc"]*out$S/(parameters["Sc"]*out$S[1])
  # out$risk <- head(c(0,out$risk),-1) # see above re index shift
  # Calculating the "interval from conversion" (see Styblo 1991 & TSRU progress report 1967)
  out$int <- parameters["Sc"]*out$S/sum(parameters["Sc"]*out$S)
  out$int <- head(c(0,out$int),-1) # see above re index shift
  # Calculating point prevalence
  out$prev <- (out$C + out$Y) / (out$E + out$Q + out$S + out$K + out$Z + out$C + out$Y)
  # Removing the first data point at time zero as we only sample after a year
  out <- out[-c(1), ]
  # Formatting the data for plotting - putting into a melted data.fame, with additional columns for the source matrix
  out <- melt(out, id.vars = c("time"))
  out$source <- source
  return(out)
}
outA <- calc(ts = times, tb = timeb, state = state, fxn = regr, parameters = paramA, source = "A")
outW <- calc(ts = times, tb = timeb, state = state, fxn = regr, parameters = paramW, source = "W")
outS <- calc(ts = times, tb = timeb, state = state, fxn = regr, parameters = paramS, source = "S")
outF <- calc(ts = times, tb = timeb, state = state, fxn = regr, parameters = paramF, source = "F")
out <- rbind(outA, outW, outS, outF)
# Plot output
theme_set(theme_bw())
##all scenarios
ggplot(out[out$variable %in% c("int"), ], aes(time, value)) + geom_point(size=2) + labs(x = "Time", y = "Incidence") + facet_grid(source ~ . , scales = "free")
##slow and slow+regression only (S+F)
out2 <- subset(out, source=='F' | source == 'S', select=time:source)  
#ggplot(out2[out2$variable %in% c("inc"), ], aes(time, value)) + geom_point(size=2) + labs(x = "Time", y = "Incidence") + facet_grid(source ~ . , scales = "free")

# Fitting parameters and state: proportion of individuals diseased after "Otime" years, using parameter set "Osource"
Osource <- "F"; Otime <- 5
Ci <- out[which(out$source == Osource & out$variable == "T" & out$time == Otime), ]
Ic <- out[which(out$source == Osource & out$variable == "inc" & out$time == Otime), ]
Ra <- out[which(out$source == Osource & out$variable == "rate" & out$time == Otime), ]
Pp <- out[which(out$source == Osource & out$variable == "prev" & out$time == Otime), ]
Rr <- out[which(out$source == Osource & out$variable == "risk" & out$time == Otime), ]
Iv <- out[which(out$source == Osource & out$variable == "int" & out$time == Otime), ]

# # Tornado plot
# range <- 0.01 ##sets range for the change in the parameters
# torn <- function(ts, tb, tt, state, fxn, parameters, source){
#   # Store data
#   data <- rbind(parameters, parameters)
#   rownames(data) <- c('+1%', '-1%')   
#   # Compare to default data set
#   out <- calc(ts = ts, tb = tb, state = state, fxn = regr, parameters = parameters, source = source)
#   def <- out[which(out$time ==tt & out$variable == "prev"),"value"]
#   for (i in 1:15){
#     # Increasing and decreasing each parameter in turn
#     parametersM = parameters; parametersL = parameters
#     parametersM[i] = parametersM[i] + range*parametersM[i]; parametersL[i] = parametersL[i] - range*parametersL[i]
#     outM <- calc(ts = ts, tb = tb, state = state, fxn = regr, parameters = parametersM, source = source); outL <- calc(ts = ts, tb = tb, state = state, fxn = regr, parameters = parametersL, source = source)
#     outM <- (outM[which(outM$time ==tt & outM$variable == "prev"),"value"] - def) / def; outL <- (outL[which(outL$time ==tt & outL$variable == "prev"),"value"] - def) / def
#     data[1, i] <- outM; data[2, i] <- outL
#   }
#   return(data)
# }
# data <- torn(ts = times, tb = timeb, tt = timet, state = state, fxn = regr, parameters = paramF, source ="F")
# # For plotting '%' on x-axis
# x <- seq(-0.01,0.01, length=10)
# ORD = order(abs(data[2,] - data[1,]))
###order black = increase in parameter, white is decrease in parameter value
# barplot(data[1,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('black'))
# barplot(data[2,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('white'), add = TRUE)
# axis(1, at=pretty(x), lab=paste0(pretty(x) * 100," %"), las=TRUE)
