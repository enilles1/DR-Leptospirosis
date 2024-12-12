######## 
# Reverse catalytic model
# NOTE
# This is the code used for the serocatalytic models presented in the manuscript. 
# This script uses age as a continuous variable. 
# However, for de-identification purposes, the dataset provided has age aggregated into a categorical variable. 
# Therefore, the provided code will need to be modified to account for this. 
########

library(tidyverse)
library(rjags)
library(loo)
library(MCMCvis)
library(binom)

source("script-serocatalytic_models/samplingFunctions.R")

## Read in cleaned lepto data (NOTE -- this will need to be modified to use aggregated data)
dat <- readRDS("data/dat_lepto.RDS")
dat_5 <- readRDS("data/dat_lepto_5.RDS")
fiveYearAgeTotals <- dat_5$total

# Extract age and seroprevalence vectors for the rjags code
# Using a bernouilli distribution within the model since have individual level data
age <- dat$age_calc
seropositive <- dat$seropositive

#################################################################################
# Defining the reverse catalytic rjags model
## Using priors informed by lepto work in Fiji for the waning value (https://journals.plos.org/plosntds/article/comments?id=10.1371/journal.pntd.0010506)
## Using a gamma distribution with shape 10 and rate 80 (or scale 0.0125)
## This distribution had a mean 0.125 for the waning parameter
## Less information regarding the FOI so using a uniform distribution of 0 to 5
#################################################################################
jcode <- "model{ 
for (i in 1:length(n.pos)){
n.pos[i] ~ dbern(seropos_est[i])
seropos_est[i] = (lambda / (lambda + delta)) * (1-exp(-(lambda+delta)*age[i])) #reverse catalytic model
}
# Priors
lambda ~ dunif(0,5) #uninformative prior
delta ~ dgamma(10,80)

}"

# Run model
mcmc.length=15000
jdat <- list(n.pos= seropositive, age=age)
jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=5000)
update(jmod)
jpos = coda.samples(jmod, c("lambda","delta"), n.iter=mcmc.length, thin = 5)

# Extract ESS and Rhat and plot MCMC chains to ensure that chains have converged and to check posteriors
MCMCsummary(jpos, round = 2)
summary(jpos)
plot(jpos)


# convert mcmc.list to a matrix and extract the overall estimates for the FOI and waning
mcmcMatrix <- as.matrix(jpos)
deltaPointEst <- mcmcMatrix[,"delta"] %>% quantile(probs=c(.5,.025,.975))
lambdaPointEst <- mcmcMatrix[,"lambda"] %>% quantile(probs=c(.5,.025,.975))

lambdaPointEst %>% round(3) %>% print
deltaPointEst %>% round(3) %>% print

#################################################################################
## Extract data for plotting
## Want to plot the model but also the model Credible intervals by age
## The FOI and waning are linked - and you cannot simply use the mean and CrI 
## Therefore:
## 1. randomly sample from the MCMC chains and the run the model for each 
## parameter selection
## 2. Then use these to obtain the CrI (age quantiles) for each age, and plot this instead
## 3. Also want to account for the sample uncertainty as well as the model uncertainty
## Therefore also take into account the sample size
#################################################################################

## Number of random samples from the MCMC chain - need to ensure it's high enough 
## to provide a good estimate
numberOfSamples = 1000
ageVector = 0:80
## Create an empty matrix where each column is age, and empty rows corresponding the number of samples
outDf <- matrix(NA,nrow=numberOfSamples, ncol = length(ageVector))


for (i in 1:numberOfSamples){
  ## Select a random number between 1 and the number of samples in the mcmc chain
  randomNumber <- floor(runif(1, min = 1, max = nrow(mcmcMatrix)))
  
  ## Extract the lambda and delta sample for that row
  ## These parameters are linked so need to extract both from the same run
  lambdaSample <- mcmcMatrix[randomNumber,2]
  deltaSample <- mcmcMatrix[randomNumber,1]
  
  ## Using these sampled parameter estimates run the model
  newRow <- (lambdaSample / (lambdaSample+deltaSample)) *(1 - exp(-ageVector*(lambdaSample+deltaSample)))
  outDf[i,] <- newRow
}

## Obtain the quantiles for each column (each one year age group). These give you
## the 95% CI
ageQuantilesModelUncertainty <- ageQuantiles(outDf)

## Extract the midpoints from the 5 year aggregated seroprevalence data
midpoints <- dat_5$midpoint

## Obtained the model uncertainty, but now want to account for size of the sample within the uncertainty estimate
## Sample from the mcmc chain as before but this time for each 5 year age group...
## Extract the binomial sampling uncertainty for this age group based on the total size of the group
randomlySampleMcmcChain <- mcmcRandomSampler(1000,mcmcMatrix,midpoints,fiveYearAgeTotals)
## Extract quantiles for the sampled results
ageQuantilesSamplingUncertainty <- ageQuantiles(randomlySampleMcmcChain)

## Add these estimates within a data frame for plotting
df_upperLower = data.frame(
  midpoint = ageVector,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-ageVector*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesModelUncertainty[,3],
  lower = ageQuantilesModelUncertainty[,2]
)

df_sampling = data.frame(
  midpoint = dat_5$midpoint,
  mean = (lambdaPointEst[1] / (lambdaPointEst[1]+deltaPointEst[1])) *(1 - exp(-dat_5$midpoint*(lambdaPointEst[1]+deltaPointEst[1]))),
  upper = ageQuantilesSamplingUncertainty[,3],
  lower = ageQuantilesSamplingUncertainty[,2]
)

## Create plot which includes model and sampling uncertainty
ggplot(df_upperLower, aes(x=midpoint, y=mean, ymin=lower, ymax=upper)) +
  geom_ribbon(alpha=0.3)+
  geom_line()+
  geom_point(data=dat_5)+
  geom_linerange(data=dat_5) +
  geom_ribbon(data=df_sampling, alpha=0.3, fill = "#457b9d")+
  scale_y_continuous(breaks=seq(0,1,by=0.1), limits = c(0,0.55))+
  xlab("Age (years)") + ylab("Proportion seropositive") +
  scale_x_continuous(breaks=seq(0,100,by=5)) +
  theme_minimal()