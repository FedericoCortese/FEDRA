# clears workspace:  
rm(list=ls()) 

# sets working directories:

setwd('C:/Users/Jeff/Documents/R/BARTmethods')

# loads necessary library
library(R2jags)
library(polspline) # helpful for Savage-Dickey analysis

#loads data (needs to be in long form; each row is a trial)
DifferenceData <- read.csv("BARTGroupDifference.csv")


# Prepare data to pass to JAGS --------------------------------------------

# Set the needed values
nParticipants = max(DifferenceData$participant) # number of participants in the experiment
maxPumps <- max(DifferenceData$burst)+ 2 # just needs to be higher than the number of pumps and participant did
totalTrials <- length(DifferenceData$participant) # total number of trials across all participants

# Specify data as vectors
burst <- DifferenceData$burst # vector of burst points on every trial
y <- DifferenceData$y # vector of participant behavior (how many times they pumped on a given trial)
participant <- DifferenceData$participant # vector of Participant ID's, indicating who is doing a given trial
z <- DifferenceData$z[1:nParticipants] # vector of group indicators. There is only 1 indicator for each participant,
                                       # so the rows don't really line up (it makes it look like participant 1's
                                       # smoker vs. nonsmoker status is flipping back and forth mid experiment).
                                       # This is an artifact of implementing the model in JASP.

data <- list("nParticipants","totalTrials","maxPumps","burst","y","participant","z") # to be passed on to JAGS

# Run model in JAGS -------------------------------------------------------

# parameters to be monitored:	
parameters <- c("delta","yPrime")

# Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# Initializes yPrime to y (ensures initial yPrime sample is possible under censoring constraint of given trial)
# Initializing sigma is not strictly necessary (but best practice for variances)
myinits <- list(list("yPrime" = y))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
samples <- jags(data, inits=myinits,parameters.to.save = parameters,
                model.file="BARTGroupDifference_1_jags.txt", n.chains=1,n.burnin = 1000, n.iter=10000, n.thin=1, DIC=T)

#### For running multiple chains in parallel. Involves different jags call and initialization ####
# When running chains in parallel, R basically initializes multiple jags calls with
# a single chain, so you build an initialization function to create new initializations
# each time R initiates the separate, parallel chain
inits.jagsParallel=function(nParticipants){
  return(list(list("yPrime" = y)))
}

#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

samples =jags.parallel(data,inits=inits.jagsParallel(nParticipants), parameters.to.save = parameters,
                       model.file="BARTGroupDifference_1_jags.txt",n.chains=3,n.burnin = 1000, n.iter=10000, n.thin=1)

#### End parallel running section ####


# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection. yPrime is the posterior predictive.

delta_samples <- samples$BUGSoutput$sims.list$delta

summary <- samples$BUGSoutput$summary

# Savage-Dickey test on r -------------------------------------------------

delta_prior = dnorm(0,0,10) # find the prior at 0

delta_post = dlogspline(0,logspline(delta_samples)) # find the posterior at 0

BF = delta_post/delta_prior

if(BF >= 1){
  BF_01 = BF
}  else{BF_10 = 1/BF}

hist(delta_samples)