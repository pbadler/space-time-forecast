
# clean up
rm(list=ls())

# set working directory
setwd("C:/Repos/space-time-forecast/genetic-diversity")

# load packages and functions
source("genetic_diversity_functions.R")


###
### 1. Simulate and describe spatial pattern in abundances at sites across a range of mean temperatures
###

sim.yrs<-2000
Tmean = seq(-5,5,0.25)   # set of mean temperatures to explore
Tstdev <- 1   # st dev of temperature

# generate one set of temperatures to use in all simulations
temperature<-matrix(NA,sim.yrs,length(Tmean))
for(iT in 1:length(Tmean)){
  temperature[,iT] = rnorm(sim.yrs,Tmean[iT],Tstdev)  # generate temperature time series
}  

# diploid model
simName <- "ss.5_noDominant"
fec.Tmu = c(-1,0,1)  # optimal temperature for genotypes AA, Aa, and aa
fec.Tsigma = rep(5,3)     # standard deviation in temperature response for genotypes AA, Aa, and aa
fec.max = c(100,100,100)  # fecundity for genotypes AA, Aa, and aa
G <- 1
seedSurv = 0.5  # survival of ungerminated seeds (same for both genotypes)
alpha = matrix(1,3,3)      # all competition coefficients (intra- and inter-) are equal for both genotypes

source("genetic_diversity_sim_spatial.R")

#CHECK: why zeros for least fit genotype's mean abundances at lowest temps but not highest?

# plot spatial figures
source("figures_spatial.R")

###
### 2. Simulate and describe temporal pattern in abundances at one site experiencing warming
###

baseline_yrs<-500  # number of yrs at baseline temperature
warming_yrs <- 100  # number of yrs with warming occurring
final_yrs <- 1000  # number of yrs at steadty-state, warmed temperature
sim_yrs <- baseline_yrs + warming_yrs + final_yrs # total number of years

Tstdev <- 1   # st dev of temperature (stationary)
baseT <- -1 # initial temperature
deltaT <- 5  # total change in temperature
Tmean <- baseT + c(rep(0,baseline_yrs),seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                   rep(deltaT,final_yrs))

# generate a temperature time series
temperature = rnorm(sim_yrs,Tmean,Tstdev)  # generate temperature time series
# plot(temperature,type="l")  # take a look

# use same parameters from diploid model above 

source("genetic_diversity_sim_temporal.R")

# FIGURES
par(mfrow=c(1,3),tcl=-0.2,mgp=c(2,0.5,0))

# plot full time series
plot(N,type="l",xlab="Time",ylab="N")

# plot temperature vs population growth during baseline period
plot(temperature_t1_baseline,r_baseline,xlab="Temperature",ylab="log Population growth rate")

# plot observed vs predicted population growth
plot(r_baseline,predict(temporal_model),xlab="Observed", ylab= "Predicted")

###
### 3. Generate forecasts for the temporal simulation (part 2)
###

source("genetic_diversity_forecast.R")


# plot time series with spatial and temporal predictions!
