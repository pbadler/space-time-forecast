
# clean up
rm(list=ls())

# set working directory
setwd("C:/repos/space-time-forecast/community-turnover")

# load packages and functions
source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")
#library("boot")
#library("latticeExtra")

###
### 1. Simulate a metacommunity: first stationary conditions, then temperature increase, then stationary again
###

# length of landscape
L_land <- 10

# temperature variables
Tmin <- 0; Tmax <- 15  # baseline temperature range
Tmean <-seq(Tmin,Tmax,length=L_land)   # baseline mean temperatures across the landscape
Tstdev <- 0.1   # st dev of temperature (interannual variation, stationary)
deltaT <- 5  # total change in mean temperature

# length of simulation
baseline_yrs<-500  # number of yrs at baseline temperature
warming_yrs <- 100  # number of yrs with warming occurring
final_yrs <- 200  # number of yrs at steady-state, warmed temperature
sim_yrs <- baseline_yrs + warming_yrs + final_yrs # total number of years

# generate one set of temperature time series to use in all simulations
temperature<-matrix(NA,sim_yrs,length(Tmean))
for(iT in 1:length(Tmean)){
  tmp_mean <- Tmean[iT] + c(rep(0,baseline_yrs),seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                          rep(deltaT,final_yrs))
  temperature[,iT] = rnorm(sim_yrs,tmp_mean,Tstdev)  # generate temperature time series
}  

# parameters to generate species pool
N <- 5      # number of species
Gmax <- 0.5   # max growth rate 
Gmin <- 0.2   # min growth rate
Lmax <- 1.5   # max sensitivity to competition
Lmin <- 0.7   # min sensitivity to competition 
Cmax <- 0.2   # max additional sensitivity to conspecific competition
Cmin <- 0.2   # min additional sensitivity to conspecific competition
d <- 10^(-3)   # fraction of offspring dispersing from home site

source("comm_turn_sim.R")


###
### 2. Simulate and describe temporal pattern in abundances at one site experiencing warming
###

baseline_yrs<-500  # number of yrs at baseline temperature
warming_yrs <- 100  # number of yrs with warming occurring
final_yrs <- 200  # number of yrs at steady-state, warmed temperature
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


###
### 3. Generate forecasts for the temporal simulation (part 2)
###

source("genetic_diversity_forecast.R")

# plot temporal figures
source("genetic_diversity_figures.R")
