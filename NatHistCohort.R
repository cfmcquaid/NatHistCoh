# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version after code-cleansing 25/05/2017
# Cohort model of TB including regression and a slow stream
### DIAGRAM ###############################################################################################################
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
### NOTES ################################################################################################################
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
### FUNCTIONS ############################################################################################################
# DISEASE DYNAMICS
regr <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # rates of change
    dE <- - 10*E
    dK <- + 10*(1 - Eq)*E + Zk*Z - (Kz + Kr)*K
    dR <- + Kr*K
    dZ <- + Kz*K + Sz*S + Qz*Q - (Zk + Zs + Zq)*Z
    dQ <- + 10*Eq*E + Sq*S + Zq*Z - (Qs + Qz)*Q
    dS <- + Qs*Q + Cs*C + Zs*Z - (Sc + Sq + Sz)*S
    dC <- + Sc*S + Yc*Y - (Cs + Cy + Cm)*C
    dY <- + Cy*C -(Yc + Ym)*Y
    dM <- + Ym*Y + Cm*C
    # return the rate of change
    list(c(dE, dK, dR, dZ, dQ, dS, dC, dY, dM))
  })
}
# BURN
burn <- function(tb, state, fxn, parameters){
  if (tb>0){
    # Calculate up until end of burn period
    out <- ode(y=state, times=seq(0,tb,by=1), func=fxn, parms=parameters)
    # Reset initial state removing indivs in Y & M
    stateb <- c(out[tb+1,2], out[tb+1,3], out[tb+1,4], out[tb+1,5], out[tb+1,6], out[tb+1,7], out[tb+1,8], Y=0, M=0)}
  else{stateb <- state}
  return(stateb)
}
# INTERVENTION
interv <- function(ti, state, fxn, parameters){
  if (ti>0){
    # Calculate up until intervention
    outi <- ode(y=state, times=seq(0,ti,by=1), func=fxn, parms=parameters)
    # Remove year 0
    outi <- outi[-1,]
    # Reset initial state moving indivs from Y to R
    statei <- c(outi[ti,2], outi[ti,3], outi[ti,4]+outi[ti,9], outi[ti,5], outi[ti,6], outi[ti,7], outi[ti,8], Y=0, outi[ti,10])}
  else{statei <- state
    # If no intervention output empty out matrix & identical state
    outi <- numeric()}
    intout <- list("out"=outi, "state"=statei)
  return(intout)
}
# SIMULATE
calc <- function(ts, tb, ti, state, fxn, parameters, source){
  # Run simulation, incorporating  burn period & intervention into disease dynamics
  # Calculate initial state after burn
  stateb <- burn(tb=tb, state=state, fxn=regr, parameters=parameters)
  # Calculate out matrix until intervention and state after intervention
  intout <- interv(ti=ti, state=stateb, fxn=regr, parameters=parameters)
  outi <- as.data.frame(intout$out)
  # Calculate after intervention
  out <- ode(y=intout$state, times=seq(0,ts-ti,by=1), func=fxn, parms=parameters)
  out <- as.data.frame(out)
  # Removing year 0
  out <- out[-1,]
  # Update times to account for intervention
  out$time <- out$time+ti
  # Merge output pre- & post-intervention
  out <- rbind(outi,out)
  # Calculate interval from conversion (see Styblo 1991 & TSRU progress report 1967)
  out$int <- parameters["Cy"]*out$C/sum(parameters["Cy"]*out$C)
  # Calculate prevalence as compartments Y + C
  out$prev <- out$Y+out$C
  # Melt data.fame for plotting
  out <- melt(out, id.vars=c("time"))
  # Include column for source of parameters behind data
  out$source <- source
  return(out)
}
# TORNADO
torn <- function(ts, tb, ti, tt, state, fxn, parameters, source){# Calculations for tornado plots
    # Storage fir results
    data <- rbind(parameters, parameters)
    rownames(data) <- c('+1%', '-1%')
    # Calculate with default data set
    out <- calc(ts=ts, tb=tb, ti=ti, state=state, fxn=regr, parameters=parameters, source=source)
    # Extract result
    def <- out[which(out$time==tt & out$variable=="int"),"value"]
    for (i in 1:16){
      # Increasing and decreasing each parameter in turn
      parametersM = parameters
      parametersL = parameters
      parametersM[i] = parametersM[i] + range*parametersM[i]
      parametersL[i] = parametersL[i] - range*parametersL[i]
      # Calculating with increased and decreased parameters
      outM <- calc(ts = ts, tb = tb, ti=ti, state = state, fxn = regr, parameters = parametersM, source = source)
      outL <- calc(ts = ts, tb = tb, ti=ti, state = state, fxn = regr, parameters = parametersL, source = source)
      # Extract results
      outM <- (outM[which(outM$time==tt & outM$variable=="int"),"value"] - def) / def
      outL <- (outL[which(outL$time==tt & outL$variable=="int"),"value"] - def) / def
      data[1, i] <- outM
      data[2, i] <- outL
    }
    return(data)
}
### CODE #################################################################################################################
library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr"); library("pryr");
# PARAMETER VALUES
  # Ax = rate from compartment A to compartment X
  # No regression  
  paramN <- c(Eq=0.10, Kr=0.10, Kz=0.01, Zk=0.00, Zq=0.01, Zs=0.00, Qz=0.00, Qs=1.50, Sz=0.00, Sq=0.00, Sc=1.00, Cs=0.00, Cy=1.00, Cm=0.10, Yc=0.00, Ym=0.50)
  # Regression
  paramR <- c(Eq=0.10, Kr=0.10, Kz=0.02, Zk=0.01, Zq=0.02, Zs=0.00, Qz=0.00, Qs=2.50, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60)
  # Proportion of change in parameters for tornado plot
  range <- 0.01
# INITIAL STATES
  state <- c(E=100000, K=0, R=0, Z=0, Q=0, S=0, C=0, Y=0, M=0)
# TIMESPAN
  # total simulation, burn period, intervention, tornado plot
  time <- list(s=10, b=1, i=2, t=5)
# CALCULATIONS
  # No regression, no intervention
  outN <- calc(ts=time$s, tb=time$b, ti=0, state=state, fxn=regr, parameters=paramN, source="N")
  # Regression, no intervention
  outR <- calc(ts=time$s, tb=time$b, ti=0, state=state, fxn=regr, parameters=paramR, source="R")
  # No regression, intervention
  outNi <- calc(ts=time$s, tb=time$b, ti=time$i, state=state, fxn=regr, parameters=paramN, source="Ni")
  # Regression, intervention
  outRi <- calc(ts=time$s, tb=time$b, ti=time$i, state=state, fxn=regr, parameters=paramR, source="Ri")
  # Tornado plot
  dat <- torn(ts=time$s, tb=time$b, ti=time$i, tt=time$t, state=state, fxn=regr, parameters=paramR, source="R")
# DATA
  # Inclusion of data for barplot
  outD <- data.frame(time=seq(1,time$s,by=1), variable=rep("D",10), value=c(.58,.24,.08,.05,.01,.01,.02,.01,0,0), source=rep("D",10))
  # Binding all results
  out <- rbind(outN, outR,outNi, outRi, outD)
### PLOTS ################################################################################################################
# DATA
plot.data <- {ggplot(out[out$variable%in%c("int")&out$source%in%c("R","N"),], aes(x=time,y=value,colour=source)) +
  # Barplot of the data
  geom_bar(data=out[out$variable %in% c("D"),], stat="identity", size=1, fill="white") +
  # Lineplot of model results
  geom_line(size=2) +
  # Line colours and legends
  scale_color_manual(values=c("black","navy","maroon"), name="Model", breaks=c("D","R","N"), labels=c("Data","Regression","No regression")) +
  # Axis labels
  labs(x="Time (years)", y="Proportion of conversions") +
  # Legend placement, size, text size
  theme_bw() +
  theme(legend.key.width=unit(2,"cm"), legend.justification=c(1,1), legend.position=c(.99,.99), legend.background=element_rect(size=0.5,linetype="solid",colour ="black"), legend.text=element_text(size=20), legend.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=25))}
# INTERVENTION
plot.interv <- {ggplot(out[out$variable%in%c("prev"),], aes(time,value,colour=source,linetype=source)) +
  # Lineplot of model results
  geom_line(size=2) +
  # Line types and legends
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed"), name="Model",breaks=c("R","N","Ri","Ni"), labels=c("Regression","No regression","Regression + intervention","No regression + intervention")) +
  # Line colours and legends
  scale_color_manual(values=c("navy","navy","maroon","maroon"),name="Model", breaks=c("R","N","Ri","Ni"), labels=c("Regression","No regression","Regression + intervention","No regression + intervention")) +
  # Axis labels
  labs(x="Time (years)", y="Prevalence (/100,000 individuals)") +
  # Axis limits
  xlim(1, 10) +
  ylim(0, 6000) +
  # Legend placement, size, text size
  theme_bw() +
  theme(legend.key.width=unit(2,"cm"), legend.justification=c(1,1), legend.position=c(.99,.99), legend.background=element_rect(size=0.5,linetype="solid",colour="black"), legend.text=element_text(size=20), legend.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=25))}
# TORNADO
{dev.control('enable')
  x <- seq(-0.01,0.01, length=10)
  ORD = order(abs(dat[2,]-dat[1,]))
  # Bar colours: black = increase in parameter value, white = decrease in parameter value
  barplot(dat[1,ORD], horiz=T, las=1, xlim=c(-0.01,0.01), xaxt='n', ylab='', beside=T, col=c('black'))
  barplot(dat[2,ORD], horiz=T, las=1, xlim=c(-0.01,0.01), xaxt='n', ylab='', beside=T, col=c('white'), add=TRUE)
  axis(1, at=pretty(x), lab=paste0(pretty(x)*100,"%"), las=TRUE)}
plot.torn = recordPlot()