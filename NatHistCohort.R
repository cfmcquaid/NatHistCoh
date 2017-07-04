# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version for model fitting 26/06/2017
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
burn <- function(tb, state, parameters){
  if (tb>0){
    # Calculate up until end of burn period
    out <- ode(y=state, times=seq(0,tb,by=1), func=regr, parms=parameters)
    # Reset initial state removing indivs in Y & M
    stateb <- c(out[tb+1,2], out[tb+1,3], out[tb+1,4], out[tb+1,5], out[tb+1,6], out[tb+1,7], out[tb+1,8], Y=0, M=0)}
  else{stateb <- state}
  return(stateb)
}
# SIMULATE
calc <- function(parameters=parameters){
  # Run simulation, incorporating  burn period into disease dynamics
  # Initial states
  state <- c(E=100000, K=0, R=0, Z=0, Q=0, S=0, C=0, Y=0, M=0)
  # total simulation time and burn period
  time <- list(s=10, b=1)
  # Calculate initial state after burn
  stateb <- burn(tb=time$b, state=state, parameters=parameters)
  # Calculate after burn
  out <- ode(y=state, times=seq(0,time$s,by=1), func=regr, parms=parameters)
  out <- as.data.frame(out)
  # Removing year 0
  out <- out[-1,]
  # Calculate interval from conversion (see Styblo 1991 & TSRU progress report 1967)
  out$int <- parameters["Cy"]*out$C/sum(parameters["Cy"]*out$C)
  # Calculate incidence
  out$inc <- parameters["Sc"]*out$S/(out$E+out$K+out$Z+out$Q+out$S)
  # Produce output for comparison with data
  data.frame(cbind(time=out$time, int=out$int, inc=out$inc))
}
# COST
regrcost <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calc(parameters)
  cost <- modCost(model=out, obs=dataINC, err="sd")
  # No standard deviation in data here as it is unknown
  return(modCost(model=out, obs=dataINT, err="sd", cost=cost))
}
# TRANSFORM
regrcost2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
  
  
  ############
  # regrcost(c(exp(lpars), Zs=0, Qz=0))}
  regrcost(c(exp(lpars), Kr=0.10, Kz=0.02, Zk=0.01, Zs=0.00, Qz=0.00, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60))}
  ############


### INPUT #################################################################################################################
library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr"); library("pryr"); library("FME");
# PARAMETER VALUES
  # Ax = rate from compartment A to compartment X
  # No regression  
  paramN <- c(Eq=0.10, Kr=0.10, Kz=0.01, Zk=0.00, Zq=0.01, Zs=0.00, Qz=0.00, Qs=1.50, Sz=0.00, Sq=0.00, Sc=1.00, Cs=0.00, Cy=1.00, Cm=0.10, Yc=0.00, Ym=0.50)
  # Regression
  paramR <- c(Eq=0.10, Kr=0.10, Kz=0.02, Zk=0.01, Zq=0.02, Zs=0.00, Qz=0.00, Qs=2.50, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60)
# DATA
  # sd gives weighting, so that the total data on eg interval since conversion = total data on incidence after 5 years
  # Data on the interval since conversion
  dataINT <- cbind(time=seq(1,10,by=1), int=c(.58,.24,.08,.05,.01,.01,.02,.01,0,0), sd=rep(10, 10))
  # Data on the incidence after 5 years
  dataINC <- cbind(time=c(0,5), inc=c(0,0.05), sd=rep(2, 2))
  dataINC <- cbind(time=c(5), inc=c(0.05), sd=rep(1, 1))
### FITTING ##############################################################################################################
  # Calculating residuals and costs
  fit <- regrcost(parameters=c(paramR))
  # Sensitivity functions (similar to tornado plot, look at L1 & L2)
  Sfun <- sensFun(regrcost, c(paramR))
  summary(Sfun)
  plot(Sfun, which=c("inc","int"), xlab="time", lwd=2)
  pairs(Sfun, which=c("inc","int"), col=c("blue","green"))
  # Collinearity
  ident <- collin(Sfun)
  ident<-ident[ !(ident$collinearity > 15), ]
  plot(ident, log="y")
  # Fix certain parameters
 
  
  ############
  # Pars <- paramR[c(1:5, 8:16)]
  Pars <- paramR[c(1, 5, 8)]
  ############
  
  
  Fit <- modFit(f=regrcost2, p=log(Pars))
  exp(coef(Fit))
  # Comparison of before and after fitting
  
  
  ############ 
  # ini <- calc(parameters = c(Pars, Zs=0, Qz=0))
  # final <- calc(parameters = c(exp(coef(Fit)), Zs=0, Qz=0))
  ini <- calc(parameters = c(Pars, Kr=0.10, Kz=0.02, Zk=0.01, Zs=0.00, Qz=0.00, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60))
  final <- calc(parameters = c(exp(coef(Fit)), Kr=0.10, Kz=0.02, Zk=0.01, Zs=0.00, Qz=0.00, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60))
  ############
  
  
  # Plot results
  par(mfrow = c(1,2))
  plot(dataINC, xlab = "time", ylab = "incidence")
    lines(ini$time, ini$inc, lty = 2)
    lines(final$time, final$inc)
  plot(dataINT, xlab = "time", ylab = "interval")
    lines(ini$time, ini$int, lty = 2)
    lines(final$time, final$int)
    legend("topright", c("data", "initial", "fitted"),lty = c(NA,2,1), pch = c(1, NA, NA))
  par(mfrow = c(1, 1))
# MCMC
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5


############ 
#MCMC <- modMCMC(f=regrcost2, p=Fit$par, niter=5000, jump=cov0, var0=var0, wvar0=0.1, updatecov=50)
MCMC <- modMCMC(f=regrcost2, p=exp(coef(Fit)), niter=5000, jump=cov0, var0=var0, wvar0=0.1, updatecov=50)
############ 


MCMC$pars <- exp(MCMC$pars)
summary(MCMC)
plot(MCMC, Full=TRUE)
pairs(MCMC, nsample=500)
