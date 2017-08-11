# CF McQuaid & RMGJ Houben
# Start 21/04/2017 ---> current version for model fitting 26/06/2017
# Cohort model of TB including regression and a slow stream
### DIAGRAM ###############################################################################################################
#                                      #
#                        M             # M - Death (due to disease)
#                       ^ ^            # W - Background mortality
#                      /   \           #
#      -> F <-> S <-> L <-> H          # F - Failed granuloma, S - Subclinical disease, L - Low intensity disease (C+S-), H - High intensity disease (C+S+, symptomatic)
#     /    \    ^      \   /           #
#    I       \  |       v v            # I - Infected, instantly a proportion Ic move to C, & 1-Ic to F
#     \        \v        D             # D - Diagnosed/treated
#      -> C <-> U                      # C - Controlled granuloma, U - Unstable granuloma,
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
###   How: play with Ic values, OR increase progression parameters
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
    dC <- + Uc*U - (Cu + Cr + w)*C
    dR <- + Cr*C
    dU <- + Cu*C + Su*S + Fu*FG - (Uc + Us + Uf + w)*U
    dF <- + Sf*S + Uf*U - (Sl + Fu + w)*FG ####Fs changed to Sl
    dS <- + Sl*FG + Ls*L + Us*U - (Sl + Sf + Su + w)*S####Fs changed to Sl
    dL <- + Sl*S + Hl*H - (Ls + Lh + Lm + Ld + w)*L
    dH <- + Lh*L -(Hl + Hm + Hd + w)*H
    dD <- + Ld*L + Hd*H
    dM <- + Hm*H + Lm*L
    dW <- + (C + U + FG + S + L + H)*w
    # return the rate of change
    list(c(dC, dR, dU, dF, dS, dL, dH, dD, dM, dW))
  })
}
# SIMULATE
#baseline
calcB <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Baseline
  #parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Regression in disease states
  #parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Clearance of infection
  #parameters <- c(parameters, Us=0.00, Fu=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Dynamic latent infection
  #parameters <- c(parameters, Us=0.00, Fu=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Regression, clearance, dynamic latency
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=1, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
    # Calculate initial state up until diagnosis
    out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
    stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
    # Calculate coprevalent cases
    out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
    stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
    # Calculate remaining time
    out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT BUT IT IS QUITE SIGNIFICANT
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
#regression
calcR <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Regression in disease states
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=1, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
  # Calculate initial state up until diagnosis
  out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
  stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
  # Calculate coprevalent cases
  out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
  stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
  # Calculate remaining time
  out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT BUT IT IS QUITE SIGNIFICANT
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
#clearance
calcC <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Clearance of infection
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=1, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
  # Calculate initial state up until diagnosis
  out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
  stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
  # Calculate coprevalent cases
  out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
  stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
  # Calculate remaining time
  out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT BUT IT IS QUITE SIGNIFICANT
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
#dynamic
calcD <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  parameters <- c(parameters, Us=0.00, Fu=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Dynamic latent infection
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=1, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
  # Calculate initial state up until diagnosis
  out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
  stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
  # Calculate coprevalent cases
  out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
  stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
  # Calculate remaining time
  out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT BUT IT IS QUITE SIGNIFICANT
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
#all
calcA <- function(parameters){
  # Fixed parameters (zero rates, non-TB mortality)
  parameters <- c(parameters, Us=0.00, Fu=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Regression, clearance, dynamic latency
  # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
  # Initial states
  state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
  # Times for different simulation periods (i.e. with different notification rates)
  time <- list(step=1, b1=1, b2=1, b3=10)
  # Run for different periods in which the notification rate changes
  # Calculate initial state up until diagnosis
  out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
  stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
  # Calculate coprevalent cases
  out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
  stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
  # Calculate remaining time
  out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
  # Compile output and remove burn period & repeated timesteps
  out <- rbind(out2, out3)
  out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2),]
  out[,"time"] <- out[,"time"] - time$b1
  out <- as.data.frame(out)
  # Calculate incidence, DO WE INCLUDE PEOPLE WHO DIE OF TB? I PRESUME NOT BUT IT IS QUITE SIGNIFICANT
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
#baseline
regrcostB <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calc(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  cost4 <- modCost(model=out, obs=dataRAG, cost=cost3, err="sd")
  # Run for each different data set
  return(cost3)
}
#regression
regrcostR <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calcR(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  cost4 <- modCost(model=out, obs=dataRAG, cost=cost3, err="sd")
  # Run for each different data set
  return(cost3)
}
#clearance
regrcostC <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calcC(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  cost4 <- modCost(model=out, obs=dataRAG, cost=cost3, err="sd")
  # Run for each different data set
  return(cost3)
}
#dynamic
regrcostD <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calcD(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  cost4 <- modCost(model=out, obs=dataRAG, cost=cost3, err="sd")
  # Run for each different data set
  return(cost3)
}
#all
regrcostA <- function(parameters){
  # A cost model comparing the output for a given parameter set to the data
  out <- calcA(parameters)
  # Run for each different data set
  cost1 <- modCost(model=out, obs=dataINC, err="sd")
  cost2 <- modCost(model=out, obs=dataINT, cost=cost1, err="sd")
  cost3 <- modCost(model=out, obs=dataCUM, cost=cost2, err="sd")
  cost4 <- modCost(model=out, obs=dataRAG, cost=cost3, err="sd")
  # Run for each different data set
  return(cost3)
}
# TRANSFORM
#baseline
regrcostB2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
  regrcostB(c(exp(lpars)))} # Baseline
#regression
regrcostR2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
  regrcostR(c(exp(lpars),Ls=0.30,Hl=0.30))} # Baseline
#clearance
regrcostC2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
regrcostC(c(exp(lpars), Cr=0.1))} # Clearance of infection
#dynamic
regrcostD2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
regrcostD(c(exp(lpars), Uc=0.5, Su=0.05, Sf=1.00))} # Dynamic latent infection
#all
regrcostA2 <- function(lpars){
  # Takes log(parameters) as input, fixes some, calculates cost
regrcost(c(exp(lpars), Cr=0.1, Uc=0.5, Su=0.05, Sf=1.00, Ls=0.30, Hl=0.30))} # Regression, clearance, dynamic latency
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
  #paramR <- c(Cu=0.5, Uf=0.05, Sl=2.00, Lh=2.00) # Baseline
  #paramR <- c(Cu=0.5, Uf=0.05, Sl=2.00, Ls=0.30, Lh=2.00, Hl=0.30) # Regression in disease states
  #paramR <- c(Cr=0.1, Cu=0.5, Uf=0.05, Sl=2.00, Lh=2.00) # Clearance of infection
  #paramR <- c(Cu=0.5, Uc=0.5, Uf=0.05, Su=0.05, Sf=1.00, Sl=2.00, Lh=2.00) # Dynamic latent infection
  #paramR <- c(Cr=0.1, Cu=0.5, Uc=0.5, Uf=0.05, Su=0.05, Sf=1.00, Sl=2.00, Ls=0.30, Lh=2.00, Hl=0.30) # Regression, clearance, dynamic latency
# DATA
  # sd gives weighting, so that the total data on eg interval since conversion = total data on incidence after 5 years
  # Data on the interval since conversion
  dataINT <- cbind(time=seq(1,10,by=1), int=c(.58,.24,.08,.05,.01,.01,.02,.01,0,0), sd=rep(10,10))
  # Data on the cumulative incidence after 5 years (DO WE INCLUDE MORTALITY FROM TB)
  dataCUM <- cbind(time=c(5), cum=c(0.05), sd=rep(1,1))
  # Data on the incidence after 5 years
  dataINC <- cbind(time=c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), inc=rep(0.00125,15), sd=rep(10,15))
  # Ragonnet Data on cumulative incidence with time
  dataRAG <- cbind(time=c(0.02739726,0.054794521,0.082191781,0.109589041,0.136986301,0.164383562,0.191780822,0.219178082,0.246575342,0.273972603,0.301369863,0.328767123,0.356164384,0.383561644,0.410958904,0.438356164,0.465753425,0.493150685,0.520547945,0.547945205,1.643835616,2.739726027,4.383561644,5.479452055), rag=c(0.011094675,0.020710059,0.027366864, 0.034763314,0.040680473,0.044378698,0.047337278,0.053254438,0.056952663,0.060650888,0.065088757,0.072485207,0.073964497,0.076923077,0.080621302,0.080621302,0.083579882,0.085798817,0.087278107,0.088017751,0.096893491,0.099852071,0.101331361,0.101331361), sd=rep(24,24))
# CALCULATION
  out <- calc(paramR)
# PARMETERS
  # p <- cbind(prog = c(paramR["Cu"], paramR["Uf"], paramR["Fs"], paramR["Sl"], paramR["Lh"], NA),
  #           reg = c(NA, paramR["Uc"], NA, paramR["Sf"], paramR["Ls"], paramR["Hl"]),
  #           other = c(paramR["Cr"], NA, NA, paramR["Su"], paramR["Lm"], paramR["Hm"]))
  # row.names(p) <- c("C", "U", "F", "S", "L", "H")
### FITTING ##############################################################################################################
# CHECKS
  # Calculating residuals and costs
  # fit <- regrcost(parameters=c(paramR))
  # # Sensitivity functions (similar to tornado plot, look at L1 & L2)
  # Sfun <- sensFun(regrcost, c(paramR))
  # # plot(Sfun, which=c("inc","int","rag","cum"), xlab="time", lwd=2)
  # # pairs(Sfun, which=c("inc","int","rag","cum"), col=c("blue","green", "red", "black"))
  # plot(Sfun, which=c("inc","int","cum"), xlab="time", lwd=2)
  # pairs(Sfun, which=c("inc","int","cum"), col=c("blue","green", "red", "black"))
  # # Collinearity identifying how many parameters are not collinear
  # ident <- collin(Sfun)
  # identY <- ident[!(ident$collinearity>15), ]
  # identX <- identY[!(identY$N<max(identY$N)), ]
  # plot(identY, log="y")
  # Fix certain parameters
  #Pars <- c(Cu=1,Uf=0.002,Sl=100,Lh=0.001)
  #Original: c(Cu=0.5, Uf=0.05, Sl=2.00, Lh=2.00)
  
  #Pars <- paramR[c(1,2,3,4)] # Baseline
  #Pars <- paramR[c(1,2,3,5)] # Regression in disease states
  #Pars <- paramR[c(2,3,4,5)] # Clearance of infection
  #Pars <- paramR[c(1,3,6,7)] # Dynamic latent infection
  #Pars <- paramR[c(2,4,7,9)] # Regression, clearance, dynamic latency
  FitB <- modFit(f=regrcostB2, p=log(c(Cu=1,Uf=0.002,Sl=100,Lh=0.001)))#,lower = log(c(0,0,0,0)), upper = log(c(10,2.3,2.5,2.4)))
  FitR <- modFit(f=regrcostR2, p=log(c(Cu=1,Uf=0.002,Sl=10,Lh=0.001)))
  FitC <- modFit(f=regrcostC2, p=log(c(Cu=1,Uf=0.002,Sl=100,Lh=0.001)))
  FitD <- modFit(f=regrcostD2, p=log(c(Cu=1,Uf=0.002,Sl=100,Lh=0.001)))
  FitA <- modFit(f=regrcostA2, p=log(c(Cu=1,Uf=0.002,Sl=100,Lh=0.001)))
  #exp(coef(FitB))
  ## Comparison of before and after fitting
  # ini <- calc(parameters = c(Pars)); final <- calc(parameters = c(exp(coef(Fit))))# Baseline
  # #ini <- calc(parameters = c(Pars, Ls=0.30, Hl=0.30)); final <- calc(parameters = c(exp(coef(Fit)), Ls=0.30, Hl=0.30))# Regression in disease states
  # #ini <- calc(parameters = c(Pars, Cr=0.1)); final <- calc(parameters = c(exp(coef(Fit)), Cr=0.1))# Clearance of infection
  # #ini <- calc(parameters = c(Pars, Uc=0.5, Su=0.05, Sf=1.00)); final <- calc(parameters = c(exp(coef(Fit)), Uc=0.5, Su=0.05, Sf=1.00))# Dynamic latent infection
  # #ini <- calc(parameters = c(Pars, Cr=0.1, Uc=0.5, Su=0.05, Sf=1.00, Ls=0.30, Hl=0.30)); final <- calc(parameters = c(exp(coef(Fit)), Cr=0.1, Uc=0.5, Su=0.05, Sf=1.00, Ls=0.30, Hl=0.30))# Regression, clearance, dynamic latency
  # ## Plot results
  #  par(mfrow = c(2,2))
  #  plot(dataINC, xlab = "time", ylab = "incidence", ylim=c(0,0.05))
  #    lines(ini$time, ini$inc, lty = 2)
  #    lines(final$time, final$inc)
  #  plot(dataINT, xlab = "time", ylab = "interval")
  #    lines(ini$time, ini$int, lty = 2)
  #    lines(final$time, final$int)
  #  plot(dataCUM, xlab = "time", ylab = "Cumulative Incidence")
  #    lines(ini$time, ini$cum, lty = 2)
  #    lines(final$time, final$cum)
  #  # plot(dataRAG, xlab = "time", ylab = "Ragonnet")
  #  #   lines(ini$time, ini$rag, lty = 2)
  #  #   lines(final$time, final$rag)
  #    legend("topright", c("data", "initial", "fitted"),lty = c(NA,2,1), pch = c(1, NA, NA))
  #  par(mfrow = c(1, 1))
  ## MCMC
  # var0 <- Fit$var_ms_unweighted
  # cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
  # MCMC <- modMCMC(f=regrcost2, p=Fit$par, niter=5000, jump=cov0, var0=var0, wvar0=0.1, updatecov=50)
  # MCMC$pars <- exp(MCMC$pars)
  # summary(MCMC)
  # plot(MCMC, Full=TRUE)
  # pairs(MCMC, nsample=500)
### INTERVENTION #################################################################################################################
   parametersB <- c(exp(coef(FitB))) # Baseline
   parametersR <- c(exp(coef(FitR)), Ls=0.30, Hl=0.30) # Regression in disease states
   parametersC <- c(exp(coef(FitC)), Cr=0.1) # Clearance of infection
   parametersD <- c(exp(coef(FitD)), Uc=0.5, Su=0.05, Sf=1.00) # Dynamic latent infection
   parametersA <- c(exp(coef(FitA)), Cr=0.1, Uc=0.5, Su=0.05, Sf=1.00, Ls=0.30, Hl=0.30) # Regression, clearance, dynamic latency
   intervB <- function(parameters){
     # Fixed parameters (zero rates, non-TB mortality)
     parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Baseline
     # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
     # Initial states
     state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
     # Times for different simulation periods (i.e. with different notification rates)
     time <- list(step=1, b1=1, b2=1, b3=10, b4=5)
     # Run for different periods in which the notification rate changes
     # Calculate initial state up until diagnosis
     out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
     stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
     # Calculate coprevalent cases
     out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
     stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
     # Calculate remaining until intervention time
     out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     stateb3 <- out3[time$b3/time$step+1, 2:ncol(out3)]
     # Calculate after intervention
     out4 <- ode(y=stateb3, times=seq(time$b1+time$b2+time$b3,time$b1+time$b2+time$b3+time$b4,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     # Compile output and remove burn period & repeated timesteps
     out <- rbind(out2, out3, out4)
     out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2, time$b3/time$step+3),]
     out[,"time"] <- out[,"time"] - time$b1
     out <- as.data.frame(out)
     # Calculate incidence
     out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
     # Calculate interval from conversion
     out$int <- (out$inc*100000)/out$D[nrow(out)]
     # Calculate cumulative incidence
     out$cum <- out$D/100000
     # Calculate cumulative incidence for Ragonnet data
     out$rag <- out$D/100000
     # Produce output for comparison with data
     return(out)
   }  
   intervR <- function(parameters){
     # Fixed parameters (zero rates, non-TB mortality)
     parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Regression in disease states
     # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
     # Initial states
     state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
     # Times for different simulation periods (i.e. with different notification rates)
     time <- list(step=1, b1=1, b2=1, b3=10, b4=5)
     # Run for different periods in which the notification rate changes
     # Calculate initial state up until diagnosis
     out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
     stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
     # Calculate coprevalent cases
     out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
     stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
     # Calculate remaining until intervention time
     out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     stateb3 <- out3[time$b3/time$step+1, 2:ncol(out3)]
     # Calculate after intervention
     out4 <- ode(y=stateb3, times=seq(time$b1+time$b2+time$b3,time$b1+time$b2+time$b3+time$b4,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     # Compile output and remove burn period & repeated timesteps
     out <- rbind(out2, out3, out4)
     out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2, time$b3/time$step+3),]
     out[,"time"] <- out[,"time"] - time$b1
     out <- as.data.frame(out)
     # Calculate incidence
     out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
     # Calculate interval from conversion
     out$int <- (out$inc*100000)/out$D[nrow(out)]
     # Calculate cumulative incidence
     out$cum <- out$D/100000
     # Calculate cumulative incidence for Ragonnet data
     out$rag <- out$D/100000
     # Produce output for comparison with data
     return(out)
   }  
   intervC <- function(parameters){
     # Fixed parameters (zero rates, non-TB mortality)
     parameters <- c(parameters, Us=0.00, Uc=0.00, Fu=0.00, Su=0.00, Sf=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Clearance of infection
     # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
     # Initial states
     state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
     # Times for different simulation periods (i.e. with different notification rates)
     time <- list(step=1, b1=1, b2=1, b3=10, b4=5)
     # Run for different periods in which the notification rate changes
     # Calculate initial state up until diagnosis
     out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
     stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
     # Calculate coprevalent cases
     out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
     stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
     # Calculate remaining until intervention time
     out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     stateb3 <- out3[time$b3/time$step+1, 2:ncol(out3)]
     # Calculate after intervention
     out4 <- ode(y=stateb3, times=seq(time$b1+time$b2+time$b3,time$b1+time$b2+time$b3+time$b4,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     # Compile output and remove burn period & repeated timesteps
     out <- rbind(out2, out3, out4)
     out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2, time$b3/time$step+3),]
     out[,"time"] <- out[,"time"] - time$b1
     out <- as.data.frame(out)
     # Calculate incidence
     out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
     # Calculate interval from conversion
     out$int <- (out$inc*100000)/out$D[nrow(out)]
     # Calculate cumulative incidence
     out$cum <- out$D/100000
     # Calculate cumulative incidence for Ragonnet data
     out$rag <- out$D/100000
     # Produce output for comparison with data
     return(out)
   }  
   intervD <- function(parameters){
     # Fixed parameters (zero rates, non-TB mortality)
     parameters <- c(parameters, Us=0.00, Fu=0.00, Hl=0.00, Ls=0.00, w=0.02, Ic=0.95, Cr=0.0, Lm=0.33, Hm=0.33) # Dynamic latent infection
     # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
     # Initial states
     state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
     # Times for different simulation periods (i.e. with different notification rates)
     time <- list(step=1, b1=1, b2=1, b3=10, b4=5)
     # Run for different periods in which the notification rate changes
     # Calculate initial state up until diagnosis
     out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
     stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
     # Calculate coprevalent cases
     out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
     stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
     # Calculate remaining until intervention time
     out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     stateb3 <- out3[time$b3/time$step+1, 2:ncol(out3)]
     # Calculate after intervention
     out4 <- ode(y=stateb3, times=seq(time$b1+time$b2+time$b3,time$b1+time$b2+time$b3+time$b4,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     # Compile output and remove burn period & repeated timesteps
     out <- rbind(out2, out3, out4)
     out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2, time$b3/time$step+3),]
     out[,"time"] <- out[,"time"] - time$b1
     out <- as.data.frame(out)
     # Calculate incidence
     out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
     # Calculate interval from conversion
     out$int <- (out$inc*100000)/out$D[nrow(out)]
     # Calculate cumulative incidence
     out$cum <- out$D/100000
     # Calculate cumulative incidence for Ragonnet data
     out$rag <- out$D/100000
     # Produce output for comparison with data
     return(out)
   }  
   intervA <- function(parameters){
     # Fixed parameters (zero rates, non-TB mortality)
     parameters <- c(parameters, Us=0.00, Fu=0.00, w=0.02, Ic=0.95, Lm=0.33, Hm=0.33) # Regression, clearance, dynamic latency
     # Run simulation, incorporating periods with different notification rates into disease dynamics, & calculate output
     # Initial states
     state <- c(C=100000*unname(parameters["Ic"]), R=0, U=0, FG=100000*(1-unname(parameters["Ic"])), S=0, L=0, H=0, D=0, M=0, W=0)
     # Times for different simulation periods (i.e. with different notification rates)
     time <- list(step=1, b1=1, b2=1, b3=10, b4=5)
     # Run for different periods in which the notification rate changes
     # Calculate initial state up until diagnosis
     out1 <- ode(y=state, times=seq(0,time$b1,by=time$step), func=regr, parms=c(parameters,Ld=0,Hd=0))
     stateb1 <- out1[time$b1/time$step+1, 2:ncol(out1)]
     # Calculate coprevalent cases
     out2 <- ode(y=stateb1, times=seq(time$b1,time$b1+time$b2,by=time$step), func=regr, parms=c(parameters,Ld=2,Hd=2))
     stateb2 <- out2[time$b2/time$step+1, 2:ncol(out2)] #DO I NEED TO SET DIAGNOSED TO ZERO HERE?
     # Calculate remaining until intervention time
     out3 <- ode(y=stateb2, times=seq(time$b1+time$b2,time$b1+time$b2+time$b3,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     stateb3 <- out3[time$b3/time$step+1, 2:ncol(out3)]
     # Calculate after intervention
     out4 <- ode(y=stateb3, times=seq(time$b1+time$b2+time$b3,time$b1+time$b2+time$b3+time$b4,by=time$step), func=regr, parms=c(parameters, Ld=1, Hd=1))
     # Compile output and remove burn period & repeated timesteps
     out <- rbind(out2, out3, out4)
     out <- out[-c(seq(1,time$b2/time$step,by=1), time$b2/time$step+2, time$b3/time$step+3),]
     out[,"time"] <- out[,"time"] - time$b1
     out <- as.data.frame(out)
     # Calculate incidence
     out$inc <- (out$D-c(0,out$D[-nrow(out)]))/100000
     # Calculate interval from conversion
     out$int <- (out$inc*100000)/out$D[nrow(out)]
     # Calculate cumulative incidence
     out$cum <- out$D/100000
     # Calculate cumulative incidence for Ragonnet data
     out$rag <- out$D/100000
     # Produce output for comparison with data
     return(out)
   }  
outB <- intervB(parametersB)
outR <- intervR(parametersR)
outC <- intervC(parametersC)
outD <- intervD(parametersD)
outA <- intervA(parametersA)
### PLOTS #################################################################################################################
  outplot <- melt(outA, id.vars=c("time"))
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
                lab = c("R - Sterilized granuloma","D - Diagnosed/Treated","M - Death (TB)","W - Death (other)"),
                leg = c("Sink size (individuals)"),
                typ = c("solid", "solid","dashed", "dashed"),
                col = c("navy","blue","maroon","red"),
                sca = 100000)}
  # SLOW
  slow <- {list(out = c("C","U"),
                lab = c("C - Controlled granuloma","U - Unstable granuloma"),
                leg = c("Slow size (individuals)"),
                typ = c("solid", "solid","dashed", "dashed"),
                col = c("navy","maroon"),
                sca = 100000)}
  # FAST
  fast <- {list(out = c("F","S","L","H"),
                lab = c("F - Failed granuloma","S - Subclinical disease","L - Low intensity disease","H - High intensity disease"),
                leg = c("Fast size (individuals)"),
                typ = c("solid", "solid","dashed", "dashed"),
                col = c("navy","blue","red","maroon"),
                sca = 2000)}
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