# Call from comm_turn_master.R. 
# Assumes many objects already created.
# Completes a couple more runs at higher dispersal rates
# and compares rate of change in model weights for the focal species.

# Uses existing species pool

###
### 1. Increase dispersal rates
###

d_baseline <- 0.01   # this should match value specified on comm_turn_master.R
d_vals <- c(d_baseline,0.05,0.1) # specify additional values of d to try
saved_weights_d <- matrix(NA, nrow=length(weight_spp),ncol=length(d_vals)) # matrix to store fitted model weights
saved_weights_d[,1] <- weight_spp  # store existing run

for(iRun in 2:length(d_vals)){
  
  d <- d_vals[iRun]
  
  # dispersal matrix (seeds disperse only to nearest neighbor locations)
  d <- d/dt
  seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
  seed_rain[seed_rain > 1] <- 0
  seed_rain <- seed_rain*(d*dt*0.5)
  diag(seed_rain) <- 1-d*dt
  seed_rain[1,1] = 1-dt*d/2  
  seed_rain[L_land,L_land] = 1-dt*d/2
  
  # Initialize community
  for(t in 2:sim_yrs){
    yr_temp <- temperature[t,]
    spxp[t,,] <- CommunityTempDis(spxp[t-1,,], tr, seed_rain, yr_temp, dt)
  } 
  spxp_mean_baseline <- apply(spxp[burnin_yrs:(burnin_yrs+baseline_yrs),,],MARGIN = c(2,3),FUN=mean)
  
  source("comm_turn_forecast.R")
  
  saved_weights_d[,iRun] <- weight_spp
  
}

###
### 2. Increase strength of selection by increasing sensitivity to interspp competition
###

# reset dispersal to baseline value
d <- d_baseline  
d <- d/dt
seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
seed_rain[seed_rain > 1] <- 0
seed_rain <- seed_rain*(d*dt*0.5)
diag(seed_rain) <- 1-d*dt
seed_rain[1,1] = 1-dt*d/2  
seed_rain[L_land,L_land] = 1-dt*d/2

c_baseline <- tr[,3]
c_factor <- c(0.5,1,1.5)  # alter interspp competition by this fraction
saved_weights_c <- matrix(NA, nrow=length(weight_spp),ncol=length(c_factor)) # matrix to store fitted model weights
saved_weights_c[,1] <- saved_weights_d[,1]  # this was the baseline run

for(iRun in 1:length(c_factor)){
  
  tr[,3] <- c_baseline*c_factor[iRun]
  
  # Initialize community
  for(t in 2:sim_yrs){
    yr_temp <- temperature[t,]
    spxp[t,,] <- CommunityTempDis(spxp[t-1,,], tr, seed_rain, yr_temp, dt)
  } 
  spxp_mean_baseline <- apply(spxp[burnin_yrs:(burnin_yrs+baseline_yrs),,],MARGIN = c(2,3),FUN=mean)
  
  source("comm_turn_forecast.R")

  saved_weights_c[,iRun] <- weight_spp
  
}
# this seems to slow the rate of change; initial communities less diverse

# write results to file
colnames(saved_weights_c) <- d_vals
write.csv(saved_weights_c,"simulations/saved_weights.csv",row.names=F)
