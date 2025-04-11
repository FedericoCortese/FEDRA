
# Install Rtools

# Install RStan
##install.packages("rstan", dependencies = TRUE)
library(rstan)

# Install hBayesDM
#install.packages("hBayesDM", dependencies=TRUE)
library(hBayesDM)

# IGT-ORL example

# Run the model with a given data.frame as df
# output <- igt_orl(
#   data = df, niter = 2000, nwarmup = 1000, nchain = 4, ncore = 4)

# Run the model with example data
# output <- igt_orl(
#   data = "example", niter = 2000, nwarmup = 1000, nchain = 4, ncore = 4)
output <- igt_orl(
  data = "example", niter = 200, nwarmup = 10, nchain = 1, ncore = 3)

head(output$rawdata)

# Visually check convergence of the sampling chains (should look like 'hairy caterpillars')
plot(output, type = "trace")

# Check Rhat values (all Rhat values should be less than or equal to 1.1)
rhat(output)

# Plot the posterior distributions of the hyper-parameters (distributions should be unimodal)
plot(output)

# Show the WAIC and LOOIC model fit estimates
printFit(output)

# Data.frame containing the summarized parameter values (as specified by indPars) for each subject
output$allIndPars

# Same but manually computed from the posterior
mean(output$parVals$K[,1])