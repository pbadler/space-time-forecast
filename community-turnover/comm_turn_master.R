
# clean up
rm(list=ls())

# set working directory
setwd("C:/repos/space-time-forecast/community-turnover")

# load packages and functions
source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")
library("boot")
#library("latticeExtra")

###
### 1. Simulate a metacommunity: first stationary conditions, then temperature increase, then stationary again
###

# length of landscape
L_land <- 20

# temperature variables
Tmin <- 0; Tmax <- 15  # baseline temperature range
Tmean <-seq(Tmin,Tmax,length=L_land)   # baseline mean temperatures across the landscape
Tstdev <- 2   # st dev of temperature (interannual variation, stationary)
deltaT <- 5  # total change in mean temperature

# length of simulation
burnin_yrs <- 2000
baseline_yrs <- 1000  # number of yrs at baseline temperature
warming_yrs <- 200  # number of yrs with warming occurring
final_yrs <- 500  # number of yrs at steady-state, warmed temperature
sim_yrs <- burnin_yrs+ baseline_yrs + warming_yrs + final_yrs # total number of years

# parameters to generate species pool
N <- 25      # number of species
Gmax <- 0.5   # max growth rate 
Gmin <- 0.2   # min growth rate
Lmax <- 1.5   # max sensitivity to competition
Lmin <- 0.7   # min sensitivity to competition 
Cmax <- 0.2   # max additional sensitivity to conspecific competition
Cmin <- 0.2   # min additional sensitivity to conspecific competition
d <- 10^(-2)   # fraction of offspring dispersing from home site

source("comm_turn_sim.R")


###
### 2. Fit models and generate forecasts
###

source("comm_turn_forecast.R")

###
### 3. Make figures
###

source("comm_turn_figures.R")




