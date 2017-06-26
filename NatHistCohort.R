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
  # INITIAL STATES
  state <- c(E=100000, K=0, R=0, Z=0, Q=0, S=0, C=0, Y=0, M=0)
  # TIMESPAN
  # total simulation, burn period
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
  return(cbind(time=out$time, int=out$int))
}
# COST
 # A cost model comparing the output for a given parameter set to the data
regrcost <- function(parameters){
  out <- calc(parameters)
  return(modCost(model=out, obs=data, err="sd"))
}
# Aotehr fitting fxn
regrcost2 <- function(lpars)
  regrcost(exp(lpars))
### CODE #################################################################################################################
library("reshape2"); library("deSolve"); library("ggplot2"); library("plyr"); library("pryr"); library("FME");
# PARAMETER VALUES
  # Ax = rate from compartment A to compartment X
  # No regression  
  paramN <- c(Eq=0.10, Kr=0.10, Kz=0.01, Zk=0.00, Zq=0.01, Zs=0.00, Qz=0.00, Qs=1.50, Sz=0.00, Sq=0.00, Sc=1.00, Cs=0.00, Cy=1.00, Cm=0.10, Yc=0.00, Ym=0.50)
  # Regression
  paramR <- c(Eq=0.10, Kr=0.10, Kz=0.02, Zk=0.01, Zq=0.02, Zs=0.00, Qz=0.00, Qs=2.50, Sz=0.01, Sq=1.00, Sc=2.00, Cs=1.00, Cy=2.00, Cm=0.10, Yc=0.10, Ym=0.60)
# DATA
  # Inclusion of data for the fitting
  data <- data.frame(time=seq(1,10,by=1), int=c(.58,.24,.08,.05,.01,.01,.02,.01,0,0), sd=rep(.1,10)) #RANDOM standard deviation chosen for now
### FITTING ##############################################################################################################
  fit <- regrcost(parameters=c(paramR))
  Sfun <- sensFun(regrcost, c(paramR))
  pairs(Sfun, which = c("int"), col = "blue")
  ident <- collin(Sfun)
  ident<-ident[ ! (ident$collinearity >15), ]
  plot(ident, log = "y")
  # fix some parameters

  Pars <- paramR[1:14] * 2
  Fit <- modFit(f = regrcost2, p = log(Pars))
  exp(coef(Fit))
  
  ini <- regr(time$s,state,Pars)
  final <- regr(time$s,state,exp(coef(Fit)))
  par(mfrow = c(1,2))
  plot(outDcomp, xlab = "time", ylab = "int")