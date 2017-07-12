# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version for model fitting 26/06/2017
# Cohort model of TB including regression and a slow stream
### DIAGRAM ###############################################################################################################
#                                      #
#                        M             # M - Death (due to disease)
#                       ^ ^            # W - Background mortality
#                      /   \           #
#      -> Q <-> S <-> C <-> Y          # Q - Failed granuloma, S - Subclinical, C - Minimal (C+S-), Y - Disseminated (C+S+, symptomatic)
#     /    \    ^      \   /           #
#    E       \  |       v v            # E - Exposed, instantly a proportion Ek move to K, & 1-Ek to Q
#     \        \v        D             # D - Diagnosed
#      -> K <-> Z                      # K - Controlled granuloma, Z- Leaky granuloma,
#         |                            #
#         v                            #
#         R                            # R - Sterilized granuloma (self-cure)
#                                      #
### NOTES ################################################################################################################
##set wd depending on who is coding/playing
setwd("C:/Users/cfmcquaid/Simulations")
#setwd("/Users/ReinHouben/Filr/ifolder/Applications (funding and jobs)/2016 - ERC starting grant/Rmodel_nathis/NatHistCoh")
# URL for fitting vignette: https://cran.r-project.org/web/packages/FME/vignettes/FME.pdf

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
    dK <- + Zk*Z - (Kz + Kr + w)*K
    dR <- + Kr*K
    dZ <- + Kz*K + Sz*S + Qz*Q - (Zk + Zs + Zq + w)*Z
    dQ <- + Sq*S + Zq*Z - (Qs + Qz + w)*Q
    dS <- + Qs*Q + Cs*C + Zs*Z - (Sc + Sq + Sz + w)*S
    dC <- + Sc*S + Yc*Y - (Cs + Cy + Cm + Cd + w)*C
    dY <- + Cy*C -(Yc + Ym + Yd + w)*Y
    dD <- + Cd*C + Yd*Y
    dM <- + Ym*Y + Cm*C
    dW <- + (K + Z + Q + S + C + Y)*w
    # return the rate of change
    list(c(dK, dR, dZ, dQ, dS, dC, dY, dD, dM, dW))
  })
}
# SIMULATE
calc <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  # Zero rates
  parameters <- c(parameters, Zs=0.00, Qz=0.00, w=0.01, Ek=0.95, Kr=0.0, Cm=0.33, Ym=0.33)
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(K=100000*unname(parameters["Ek"]), R=0, Z=0, Q=100000*(1-unname(parameters["Ek"])), S=0, C=0, Y=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=0.01, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
    # Calculate initial state up until diagnosis
    out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Cd=0,Yd=0))
    stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
    # Calculate coprevalent cases
    out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Cd=2,Yd=2))
    stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)]
    # Calculate remaining time
    out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Cd=1, Yd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT
  out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
  # Calculate interval from conversion
  out$int <- (out$inc*100000)/out$D[nrow(out)]
  # Calculate cumulative incidence (DO WE DIVIDE BY 100000 OR DO WE EXCLUDE COPREVALENT CASES)
  out$cum <- out$D/100000
  # Calculate cumulative incidence for Ragonnet data
  out$rag <- out$D/100000
  # Produce output for comparison with data
  return(out)
}
# COST
regrcost <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calc(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  # Run for each different data set
  return(modCost(model=out, obs=dataRAG, cost=cost3, err="sd"))
}
# TRANSFORM
regrcost2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
  regrcost(c(exp(lpars), Kz=0.02, Qs=2.50, Sz=0.01, Sq=1.00, Cs=1.00, Cy=2.00, Yc=0.1))}
# TORNADO PLOT
torn <- function(time, variable, range, parameters){# Calculations for tornado plots
  # Storage for results
  res <- rbind(parameters, parameters)
  rownames(res) <- c('+1%', '-1%')
  # Calculate with default data set
  out <- calc(parameters=parameters)
  # Extract result
  def <- out[time, variable]
  for (i in 1:length(parameters)){
    # Increasing and decreasing each parameter in turn
    parametersM = parameters
    parametersL = parameters
    parametersM[i] = parametersM[i] + range*parametersM[i]
    parametersL[i] = parametersL[i] - range*parametersL[i]
    # Calculating with increased and decreased parameters
    outM <- calc(parameters = parametersM)
    outL <- calc(parameters = parametersL)
    # Extract results
    outM <- (outM[time, variable] - def) / def
    outL <- (outL[time, variable] - def) / def
    res[1, i] <- outM
    res[2, i] <- outL
  }
  return(res)
}
### INPUT #################################################################################################################
library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr"); library("pryr"); library("FME");
# PARAMETER VALUES
  # Ax = rate from compartment A to compartment X
  # No regression  
  paramN <- c(Kz=0.01, Zk=0.00, Zq=0.01, Qs=1.50, Sz=0.00, Sq=0.00, Sc=1.00, Cs=0.00, Cy=1.00, Yc=0.00)
  # Regression
  paramR <- c(Kz=0.02, Zk=0.01, Zq=0.02, Qs=2.50, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Yc=0.10)
# DATA
  # sd gives weighting, so that the total data on eg interval since conversion = total data on incidence after 5 years
  # Data on the interval since conversion
  dataINT <- cbind(time=seq(1,10,by=1), int=c(.58,.24,.08,.05,.01,.01,.02,.01,0,0), sd=rep(10,10))
  # Data on the cumulative incidence after 5 years (DO WE INCLUDE MORTALITY FROM TB)
  dataCUM <- cbind(time=c(5), cum=c(0.05), sd=rep(1,1))
  # Data on the incidence after 5 years
  dataINC <- cbind(time=c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), inc=rep(0.00125,15), sd=rep(15,15))
  # Ragonnet Data on incidence with time
  dataRAG <- cbind(time=c(0.02739726,0.054794521,0.082191781,0.109589041,0.136986301,0.164383562,0.191780822,0.219178082,0.246575342,0.273972603,0.301369863,0.328767123,0.356164384,0.383561644,0.410958904,0.438356164,0.465753425,0.493150685,0.520547945,0.547945205,1.643835616,2.739726027,4.383561644,5.479452055), rag=c(0.011094675,0.020710059,0.027366864, 0.034763314,0.040680473,0.044378698,0.047337278,0.053254438,0.056952663,0.060650888,0.065088757,0.072485207,0.073964497,0.076923077,0.080621302,0.080621302,0.083579882,0.085798817,0.087278107,0.088017751,0.096893491,0.099852071,0.101331361,0.101331361), sd=rep(24,24))
# CALCULATION
  out <- calc(paramR)
# PARMETERS
  p <- cbind(prog = c(paramR["Kz"], paramR["Zq"], paramR["Qs"], paramR["Sc"], paramR["Cy"], NA),
            reg = c(NA, paramR["Zk"], NA, paramR["Sq"], paramR["Cs"], paramR["Yc"]),
            other = c(paramR["Kr"], NA, NA, paramR["Sz"], paramR["Cm"], paramR["Ym"]))
  row.names(p) <- c("K", "Z", "Q", "S", "C", "Y")
### PLOTS #################################################################################################################
outplot <- melt(out, id.vars=c("time"))
# Choose outputs to compare, labels, line colours & types, y-legend, scale
# INCIDENCE
inc <- {list(out = c("inc"),
  lab = c("Incidence"),
  leg = c("Incidence (proportion)"),
  typ = c("solid"),
  col = c("navy"),
  sca = 0.1)}
# INTERVAL
int <- {list(out = c("int"),
  lab = c("Interval"),
  leg = c("Proportion"),
  typ = c("solid"),
  col = c("navy"),
  sca = 1)}
# SINK
sink <- {list(out = c("R","D","M","W"),
  lab = c("R - Sterilized granuloma","D - Diagnosed","M - Death (TB)","W - Death (other)"),
  leg = c("Sink size (individuals)"),
  typ = c("solid", "solid","dashed", "dashed"),
  col = c("navy","blue","maroon","red"),
  sca = 100000)}
# SLOW
slow <- {list(out = c("K","Z"),
  lab = c("K - Controlled granuloma","Z - Leaky granuloma"),
  leg = c("Slow size (individuals)"),
  typ = c("solid", "solid","dashed", "dashed"),
  col = c("navy","maroon"),
  sca = 100000)}
# FAST
fast <- {list(out = c("Q","S","C","Y"),
  lab = c("Q - Failed granuloma","S - Subclinical","C - Minimal","Y - Disseminated"),
  leg = c("Fast size (individuals)"),
  typ = c("solid", "solid","dashed", "dashed"),
  col = c("navy","blue","red","maroon"),
  sca = 1000)}
# PLOT
comp <- sink
plot.out <- {ggplot(outplot[outplot$variable%in%comp$out,], aes(time,value,colour=variable,linetype=variable)) +
  # Lineplot of model results
  geom_line(size=2) +
  # Line types and legends
  scale_linetype_manual(values=comp$typ, name="Model", breaks=comp$out, labels=comp$lab) +
  # Line colours and legends
  scale_color_manual(values=comp$col, name="Model", breaks=comp$out, labels=comp$lab) +
  # Axis labels
  labs(x="Time (years)", y=comp$leg) +
  # Axis limits
  xlim(1, max(outplot$time)) +
  ylim(0, comp$sca) +
  # Legend placement, size, text size
  theme_bw() +
  theme(legend.key.width=unit(2,"cm"), legend.justification=c(1,1), legend.position=c(.99,.99), legend.background=element_rect(size=0.5,linetype="solid",colour="black"), legend.text=element_text(size=20), legend.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=25))}
plot.out
### TORNADO  ###############################################################################################################
dat <- torn(time=5, variable="int", range=0.01, parameters=paramR)
{dev.control('enable')
  x <- seq(-0.01,0.01, length=10)
  ORD = order(abs(dat[2,]-dat[1,]))
  # Bar colours: black = increase in parameter value, white = decrease in parameter value
  barplot(dat[1,ORD], horiz=T, las=1, xlim=c(-0.02,0.02), xaxt='n', ylab='', beside=T, col=c('black'))
  barplot(dat[2,ORD], horiz=T, las=1, xlim=c(-0.02,0.02), xaxt='n', ylab='', beside=T, col=c('white'), add=TRUE)
  axis(1, at=pretty(x), lab=paste0(pretty(x)*100,"%"), las=TRUE)}
plot.torn = recordPlot()
### FITTING ##############################################################################################################
# CHECKS
  # Calculating residuals and costs
  fit <- regrcost(parameters=c(paramR))
  # Sensitivity functions (similar to tornado plot, look at L1 & L2)
  Sfun <- sensFun(regrcost, c(paramR))
  summary(Sfun)
  # plot(Sfun, which=c("inc","int"), xlab="time", lwd=2)
  # pairs(Sfun, which=c("inc","int"), col=c("blue","green"))
  # Collinearity identifying how many parameters are not collinear
  ident <- collin(Sfun)
  identY <- ident[!(ident$collinearity>15), ]
  identX <- identY[!(identY$N<max(identY$N)), ]
  # plot(identY, log="y")
  # Fix certain parameters
  Pars <- paramR[c(2, 3, 7)]
  Fit <- modFit(f=regrcost2, p=log(Pars))
  exp(coef(Fit))
  ## Comparison of before and after fitting
  # ini <- calc(parameters = c(Pars, Zk=0.01, Zq=0.02, Sz=0.01, Sc=2.00, Cs=1.00, Cy=2.00, Yc=0.10))
  # final <- calc(parameters = c(exp(coef(Fit)), Zk=0.01, Zq=0.02, Sz=0.01, Sc=2.00, Cs=1.00, Cy=2.00, Yc=0.10))
  ## Plot results
  # par(mfrow = c(1,2))
  # plot(dataINC, xlab = "time", ylab = "incidence")
  #   lines(ini$time, ini$inc, lty = 2)
  #   lines(final$time, final$inc)
  # plot(dataINT, xlab = "time", ylab = "interval")
  #   lines(ini$time, ini$int, lty = 2)
  #   lines(final$time, final$int)
  #   legend("topright", c("data", "initial", "fitted"),lty = c(NA,2,1), pch = c(1, NA, NA))
  # par(mfrow = c(1, 1))
  ## MCMC
  var0 <- Fit$var_ms_unweighted
  cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
  MCMC <- modMCMC(f=regrcost2, p=Fit$par, niter=5000, jump=cov0, var0=var0, wvar0=0.1, updatecov=50)
  MCMC$pars <- exp(MCMC$pars)
  summary(MCMC)
  plot(MCMC, Full=TRUE)
  pairs(MCMC, nsample=500)