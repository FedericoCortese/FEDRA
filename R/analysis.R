# Load and clean data -----------------------------------------------------

load("FEDRA_varsxclust_all.RData")
load("FEDRA_varsxclust_k12.RData")
load("FEDRA_varsxclust_k7.RData")

load("FEDRA_data_iowa.RData")
load("FEDRA_data_bart.RData")

# In ogni file c'è una variabile IDK che rappresenta il soggetto e una variabile GRP che rappresenta il gruppo. 
# I soggetti con GRP = EXC vanno esclusi.

exc=data[which(data$GRP=="EXC"),1]

data_cleaned=data[-which(data$IDK%in%exc),]
data_k7_cleaned=data_k7[-which(data_k7$IDK%in%exc),]
data_k12_cleaned=data_k12[-which(data_k12$IDK%in%exc),]

data_iowa_cleaned=data_iowa[-which(data_iowa$IDK%in%exc),]
data_bart_cleaned=data_bart[-which(data_bart$IDK%in%exc),]

# save(data_cleaned,data_k7_cleaned,data_k12_cleaned,
#      data_iowa_cleaned,data_bart_cleaned,
#      file="FEDRA_data_cleaned.RData")


# Load cleaned data -------------------------------------------------------

load("FEDRA_data_cleaned.RData")

# IGT ---------------------------------------------------------------------

# Input data format
# trial choice gain loss subjID
# 1     1      3   50    0   1001
# 2     2      2  100    0   1001
# ...  ...    ... ...   ...  ...

# Structure: all int

head(data_iowa_cleaned)
str(data_iowa_cleaned)

data_iowa_cleaned$trial=as.integer(data_iowa_cleaned$trial)
data_iowa_cleaned$choice=as.integer(
  dplyr::recode(data_iowa_cleaned$deck, "A"=1, "B"=2, "C"=3, "D"=4))
data_iowa_cleaned$win=as.integer(data_iowa_cleaned$win)
data_iowa_cleaned$lose=as.integer(data_iowa_cleaned$lose)

data_iowa_cleaned$subjID=as.integer(gsub("#", "", data_iowa_cleaned$IDK))

#data_iowa_final=data_iowa_cleaned[,c("trial", "choice", "win", "lose", "subjID")]
data_iowa_final=data_iowa_cleaned[,c("subjID","choice", "win", "lose")]
colnames(data_iowa_final)=c("subjID","choice","gain","loss")

head(data_iowa_final)
tail(data_iowa_final)
str(data_iowa_final)

write.table(data_iowa_final, file = "data_iowa_final.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#install.packages("hBayesDM", dependencies=TRUE)

# output_iowa <- hBayesDM::igt_orl(
#   data = "data_iowa_final.txt", niter = 400, nwarmup = 100, nchain = 1
#   )

library(hBayesDM)
library(rstan)

output_iowa <- hBayesDM::igt_orl(
  data = "data_iowa_final.txt", niter = 4000, nwarmup = 1000, nchain = 3, ncore = parallel::detectCores()-1)
#save(output_iowa, file="output_iowa.RData")

# Number of subjects
length(unique(data_iowa_final$subjID))

# Trace plot for K for the first subject
plot(output_iowa$parVals$K[,1],type='l')

# Check Rhat values (all Rhat values should close to 1)
Rhat_iowa=rhat(output_iowa)
head(Rhat_iowa)
hist(Rhat_iowa$Rhat)
summary(Rhat_iowa$Rhat)

# Means of the posteriors
data_iowa_posteriors=output_iowa$allIndPars
data_iowa_posteriors$IDK=unique(data_iowa_cleaned$IDK)
save(data_iowa_posteriors, file="data_iowa_posteriors.RData")

# BART --------------------------------------------------------------------

# Input data format
# participant y burst
#  1          1     1
#  1          1    13
#  1          2     7
# ...        ...   ... 

# Structure: all integers

head(data_bart_cleaned)

length(unique(data_bart_cleaned$IDK))

data_bart_cleaned$participant=as.integer(gsub("#", "", data_bart_cleaned$IDK))
data_bart_cleaned$y=as.integer(data_bart_cleaned$pumps)
data_bart_cleaned$burst=as.integer(data_bart_cleaned$life)

data_bart_final=data_bart_cleaned[,c("participant", "y", "burst")]

head(data_bart_final)
tail(data_bart_final)

# Set the needed values
nParticipants = length(unique(data_bart_final$participant)) # number of participants in the experiment
maxPumps <- max(data_bart_final$burst)+ 2 # just needs to be higher than the number of pumps and participant did
totalTrials <- length(data_bart_final$participant) # total number of trials across all participants

burst <- data_bart_final$burst # vector of burst points on every trial
y <- data_bart_final$y # vector of participant behavior (how many times they pumped on a given trial)
participant <- data_bart_final$participant # vector of Participant ID's, indicating who is doing a given trial

data <- list("nParticipants","totalTrials","maxPumps","burst","y","participant") # to be passed on to JAGS

# parameters to be monitored:	
parameters <- c("beta", "rho","yPrime")

# Fixed initializations (Have to have a list in a list for jags, for more than one chain, add another sub list)
# Initializes yPrime to y (ensures initial yPrime sample is possible under censoring constraint of given trial)
# Initializing sigma is not strictly necessary (but best practice for variances)
myinits <- list(list("yPrime" = y,"beta" = runif(nParticipants,.5,1)))

# The following command calls JAGS with specific options.
#This is for running 1 chain (use code below for faster multiple chains)
# samples <- jags(data, inits=myinits,parameters.to.save = parameters,
#                 model.file="BARTMeasurement_1_jags.txt", 
#                 n.chains=1,
#                 n.burnin = 1000, 
#                 n.iter=10000, n.thin=1, DIC=T)

#### For running multiple chains in parallel. Involves different jags call and initialization ####
# When running chains in parallel, R basically initializes multiple jags calls with
# a single chain, so you build an initialization function to create new initializations
# each time R initiates the separate, parallel chain
inits.jagsParallel=function(nParticipants){
  return(list(list("yPrime" = y,"beta" = runif(nParticipants,.5,1))))
}

#Now use jags.parallel to run multiple chains much quicker, adjust chains in n.chains

output_bart =R2jags::jags.parallel(data,inits=inits.jagsParallel(nParticipants), 
                       parameters.to.save = parameters,
                       model.file="BARTMeasurement_1_jags.txt",n.chains=3,
                       n.burnin = 1000, n.iter=10000, 
                       n.thin=1)

#save(output_bart, file="output_bart.RData")

rho_samples <- output_bart$BUGSoutput$sims.list$rho
# x11()
# matplot(output_bart$BUGSoutput$sims.list$rho,type='l')

beta_samples <- output_bart$BUGSoutput$sims.list$beta
x11()
matplot(output_bart$BUGSoutput$sims.list$beta,type='l')


summary <- output_bart$BUGSoutput$summary
summary

av_betas=summary[1:nParticipants,1]

dim(output_bart$BUGSoutput$sims.matrix)
