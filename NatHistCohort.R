# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version after revision with Hanif 18/05/2017
# Cohort model of TB including regression and a slow stream
########################################
#                                      #
#                        M             # M - Death (due to disease)
#                       ^ ^            #
#                      /   \           #
#      -> Q <-> S <-> C <-> Y          # Q - Failed granuloma, S - Subclinical, C - Minimal (C+S-), Y - Disseminated (C+S+, symptomatic)
#     /    \    ^                      #
#    E       \  |                      # E - Exposed
#     \        \v                      #
#      -> K <-> Z                      # K - Controlled granuloma, Z- Leaky granuloma,
#         |                            #
#         v                            #
#         R                            # R - Sterilized granuloma (self-cure)
#                                      #
########################################
##set wd depending on who is coding/playing
setwd("C:/Users/cfmcquaid/Simulations")
#setwd("/Users/ReinHouben/Filr/ifolder/Applications (funding and jobs)/2016 - ERC starting grant/Rmodel_nathis/NatHistCoh")

### URL for how to upload changes: http://r-bio.github.io/intro-git-rstudio/
### to pull in Finn's version - type in shell: git pull upstream master
### to push my version - create pull request to Finn:  
#       commit (remember to tick box), then push
#       goto Github, pull requests, create pull request, if merge can be done on auto, just complete
#       Alternatively Finn complete pull request (he does not get email notification)
# Finn to just click push. Job done

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


# No regression: 1/Cy + 1/Ym = 3 years, Ym = 0.5 (0.3 death 0.2 self-cure), therefore Cy = 1. NOTE: 70% die, 30% self-cure
# Regression: 1/(Yc + Ym) = 0.5 years, Ym = 0.5 (as above), therefore Yc = 1.5.
# Regression: 1/(Cy + Cs) = 0.5 years. Assume ratio 2:1 for Cy:Cs, therefore Cs = 2/3
paramS <- c(Eq=0.10, Qs=1.50, Qz=0.00, Kz=0.01, Kr=0.10, Sc=1.00, Sq=0.00, Sz=0.00, Zk=0.00, Zs=0.00, Zq=0.01, Cy=1.00, Cs=0.00, Cm=0.10, Ym=0.50, Yc=0.00)
paramF <- c(Eq=0.10, Qs=2.50, Qz=0.00, Kz=0.02, Kr=0.10, Sc=2.00, Sq=1.00, Sz=0.01, Zk=0.01, Zs=0.00, Zq=0.02, Cy=2.00, Cs=1.00, Cm=0.10, Ym=0.60, Yc=0.10)
# Initial states of compartments
state <- c(E=100000, Q=0, K=0, S=0, Z=0, C=0, Y=0, R=0, M=0)
# Timespan for simulation
times <- seq(0, 10, by = 1)
# Timespan for burn
timeb <- 1
# Timespan for intervention
timei <- 2
# Timespan for tornado
timet <- 5
# ODE function
regr <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
   # rates of change
   dE <- - 10*E
   dQ <- + 10*Eq*E + Sq*S + Zq*Z - (Qs + Qz)*Q
   dK <- + 10*(1 - Eq)*E + Zk*Z - (Kz + Kr)*K
   dS <- + Qs*Q + Cs*C + Zs*Z - (Sc + Sq + Sz)*S
   dZ <- + Kz*K + Sz*S + Qz*Q - (Zk + Zs + Zq)*Z
   dC <- + Sc*S + Yc*Y - (Cs + Cy + Cm)*C
   dY <- + Cy*C -(Yc + Ym)*Y
   dR <- + Kr*K
   dM <- + Ym*Y + Cm*C
   # return the rate of change
   list(c(dE, dQ, dK, dS, dZ, dC, dY, dR, dM))
  })
}
# Burn off the first tb years to establish the proportion of individuals in each compartment, removing diseased
burn <- function(tb, state, fxn, parameters){
  if (tb>0){
  out <- ode(y = state, times = seq(0, tb, by = 1), func = fxn, parms = parameters)
  stateb <- c(out[tb+1, 2], out[tb+1, 3], out[tb+1, 4], out[tb+1, 5], out[tb+1, 6], out[tb+1, 7], Y=0, out[tb+1, 9], M=0)}
  else{stateb <- state}
  return(stateb)
}
# An intervention after ti years moves individuals from Y to R
interv <- function(ti, state, fxn, parameters){
  if (ti>0){
    outi <- ode(y = state, times = seq(0, ti, by = 1), func = fxn, parms = parameters)
    outi <- outi[-1,]
    statei <- c(outi[ti, 2], outi[ti, 3], outi[ti, 4], outi[ti, 5], outi[ti, 6], outi[ti, 7], Y=0, outi[ti, 9]+outi[ti, 8], outi[ti, 10])}
  else{statei <- state
  outi <- numeric()}
  intout <- list("out" = outi,"state" = statei)
  return(intout)
}
# Calculation function
calc <- function(ts, tb, ti, state, fxn, parameters, source){
  stateb <- burn(tb = tb, state = state, fxn = regr, parameters = parameters)
  intout <- interv(ti=ti, state = stateb, fxn = regr, parameters = parameters)
  outi <- as.data.frame(intout$out)
  out <- ode(y = intout$state, times = ts, func = fxn, parms = parameters)
  out <- as.data.frame(out)
  # Removing the first timestep, time zero, and all timesteps after "times"
  out <- out[-1,]
  out$time <- out$time+ti
  out <- rbind(outi,out)
  out <- out[1:10,]
  # Calculating the "interval from conversion" (see Styblo 1991 & TSRU progress report 1967)
  out$int <- parameters["Cy"]*out$C/sum(parameters["Cy"]*out$C)
  out$prev <- out$Y+out$C
  # Formatting the data for plotting - putting into a melted data.fame, with additional columns for the source matrix
  out <- melt(out, id.vars = c("time"))
  out$source <- source
  return(out)
}
# No intervention
outS <- calc(ts = times, tb = timeb, ti=0, state = state, fxn = regr, parameters = paramS, source = "S")
outF <- calc(ts = times, tb = timeb, ti=0, state = state, fxn = regr, parameters = paramF, source = "F")
# Intervention
outSi <- calc(ts = times, tb = timeb, ti=timei, state = state, fxn = regr, parameters = paramS, source = "Si")
outFi <- calc(ts = times, tb = timeb, ti=timei, state = state, fxn = regr, parameters = paramF, source = "Fi")
# Including the data for barplot
time <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10); variable <- c("D", "D", "D", "D", "D", "D", "D", "D", "D", "D"); value <- c(0.58, 0.24, 0.08, 0.05, 0.01, 0.01, 0.02, 0.01, 0, 0); source <- c("D", "D", "D", "D", "D", "D", "D", "D", "D", "D")
outD <- data.frame(time, variable, value, source)
out <- rbind(outS, outF,outSi, outFi, outD)
# Plot output
ggplot(out[out$variable %in% c("prev"), ], aes(time, value, colour=source, linetype=source)) + geom_line(size=2) + scale_linetype_manual(values=c("solid", "dashed","solid", "dashed"),name="Model",breaks=c("F","S","Fi","Si"),labels=c("Regression", "No regression", "Regression + intervention", "No regression + intervention")) + scale_color_manual(values=c("navy","navy","maroon","maroon"),name="Model",breaks=c("F","S","Fi","Si"),labels=c("Regression", "No regression", "Regression + intervention", "No regression + intervention")) + labs(x = "Time (years)", y = "Prevalence (/100,000 individuals)") + xlim(1,10)+ ylim(0,6000)+theme_bw()+ theme(legend.key.width=unit(2,"cm"),legend.justification=c(1,1),legend.position=c(.99,.99),legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),legend.text=element_text(size=20),legend.title=element_text(size=20),axis.text = element_text(size = 15),axis.title = element_text(size = 25))#+ facet_wrap(~ variable , scales = "free", nrow = 2)
ggplot(out[out$variable %in% c("int")&out$source %in% c("F","S"), ], aes(x = time, y = value, colour = source)) + geom_bar(data=out[out$variable %in% c("D"), ], stat="identity",size=1,fill="white") + geom_line(size=2) + scale_color_manual(values=c("black","navy","maroon"),name="Model",breaks=c("D","F","S"),labels=c("Data","Regression", "No regression")) + labs(x = "Time (years)", y = "Proportion of conversions") + theme_bw() + theme(legend.key.width=unit(2,"cm"),legend.justification=c(1,1),legend.position=c(.99,.99),legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),legend.text=element_text(size=20),legend.title=element_text(size=20),axis.text = element_text(size = 15),axis.title = element_text(size = 25))
# Fitting parameters and state: proportion of individuals diseased after "Otime" years, using parameter set "Osource"
Osource <- "F"; Otime <- 5
Iv <- out[which(out$source == Osource & out$variable == "int" & out$time == Otime), ]


# # # Tornado plot
# range <- 0.01 ##sets range for the change in the parameters
# torn <- function(ts, tb, ti, tt, state, fxn, parameters, source){
#   # Store data
#   data <- rbind(parameters, parameters)
#   rownames(data) <- c('+1%', '-1%')
#   # Compare to default data set
#   out <- calc(ts = ts, tb = tb, ti=ti, state = state, fxn = regr, parameters = parameters, source = source)
#   def <- out[which(out$time ==tt & out$variable == "int"),"value"]
#   for (i in 1:16){
#     # Increasing and decreasing each parameter in turn
#     parametersM = parameters; parametersL = parameters
#     parametersM[i] = parametersM[i] + range*parametersM[i]; parametersL[i] = parametersL[i] - range*parametersL[i]
#     outM <- calc(ts = ts, tb = tb, ti=ti, state = state, fxn = regr, parameters = parametersM, source = source); outL <- calc(ts = ts, tb = tb, ti=ti, state = state, fxn = regr, parameters = parametersL, source = source)
#     outM <- (outM[which(outM$time ==tt & outM$variable == "int"),"value"] - def) / def; outL <- (outL[which(outL$time ==tt & outL$variable == "int"),"value"] - def) / def
#     data[1, i] <- outM; data[2, i] <- outL
#   }
#   return(data)
# }
# data <- torn(ts = times, tb = timeb, ti=timei, tt = timet, state = state, fxn = regr, parameters = paramF, source ="F")
# # For plotting '%' on x-axis
# x <- seq(-0.01,0.01, length=10)
# ORD = order(abs(data[2,] - data[1,]))
# ##order black = increase in parameter, white is decrease in parameter value
# barplot(data[1,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('black'))
# barplot(data[2,ORD], horiz = T, las=1, xlim = c(-0.01,0.01), xaxt='n', ylab = '', beside=T, col=c('white'), add = TRUE)
# axis(1, at=pretty(x), lab=paste0(pretty(x) * 100," %"), las=TRUE)
