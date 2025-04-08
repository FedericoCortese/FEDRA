# clears workspace:  
rm(list=ls()) 

# sets working directories:

# setwd('C:/Users/Jeff/Documents/R/BARTmethods')

# loads necessary library
library(R2jags)

#loads data (needs to be in long form; each row is a trial)
DemoData <- read.csv("balloonRiskGuanEtAl.csv")


# Prepare data to pass to JAGS --------------------------------------------

# Set the needed values
nParticipants = max(DemoData$participant) # number of participants in the experiment
maxPumps <- max(DemoData$burst)+ 2 # just needs to be higher than the number of pumps and participant did
totalTrials <- length(DemoData$participant) # total number of trials across all participants

# Specify data as vectors
burst <- DemoData$burst # vector of burst points on every trial
y <- DemoData$y # vector of participant behavior (how many times they pumped on a given trial)
participant <- DemoData$participant # vector of Participant ID's, indicating who is doing a given trial

data <- list("nParticipants","totalTrials","maxPumps","burst","y","participant") # to be passed on to JAGS

# Run model in JAGS -------------------------------------------------------

# parameters to be monitored:	
parameters <- c("beta", "rho","yPrime")

# Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# Initializes yPrime to y (ensures initial yPrime sample is possible under censoring constraint of given trial)
# Initializing sigma is not strictly necessary (but best practice for variances)
myinits <- list(list("yPrime" = y,"beta" = runif(nParticipants,.5,1)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, inits=myinits,parameters.to.save = parameters,
                model.file="BARTMeasurement_1_jags.txt", n.chains=1,n.burnin = 1000, n.iter=10000, n.thin=1, DIC=T)

#### For running multiple chains in parallel. Involves different jags call and initialization ####
# When running chains in parallel, R basically initializes multiple jags calls with
# a single chain, so you build an initialization function to create new initializations
# each time R initiates the separate, parallel chain
inits.jagsParallel=function(nParticipants){
  return(list(list("yPrime" = y,"beta" = runif(nParticipants,.5,1))))
}

#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

samples =jags.parallel(data,inits=inits.jagsParallel(nParticipants), parameters.to.save = parameters,
                       model.file="BARTMeasurement_1_jags.txt",n.chains=3,n.burnin = 1000, n.iter=10000, n.thin=1)

#### End parallel running section ####


# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection. yPrime is the posterior predictive.

rho_samples <- samples$BUGSoutput$sims.list$rho
beta_samples <- samples$BUGSoutput$sims.list$beta
  
summary <- samples$BUGSoutput$summary
