
# Call from genetic_diversity_master.R, with all objects in memory

# See how rate of shift in models weights responds to
# an increased temperature niche of each genotype.

fec_Tsigma_baseline <- fec_Tsigma # baseline temperature niche widths
delta_fec_sigma <- c(0.5,1,2)  
saved_weights <- matrix(NA,nrow=length(weight),ncol=length(delta_fec_sigma))
saved_weights[,1] <- weight

for(iRun in 1:length(delta_fec_sigma)){
  
  fec_Tsigma <- fec_Tsigma_baseline*delta_fec_sigma[iRun]
  
  # simulate baseline populations across temperature gradient
  sim_yrs<-2000
  Tmean = seq(-5,5,0.25)   # set of mean temperatures to explore
  Tstdev <- 1   # st dev of temperature
  temperature<-matrix(NA,sim_yrs,length(Tmean))
  for(iT in 1:length(Tmean)){
    temperature[,iT] = rnorm(sim_yrs,Tmean[iT],Tstdev)  # generate temperature time series
  } 
  source("genetic_diversity_sim_spatial.R")
  
  # simulate long-term run for focal population
  baseline_yrs<-500  # number of yrs at baseline temperature
  warming_yrs <- 100  # number of yrs with warming occurring
  final_yrs <- 200  # number of yrs at steadty-state, warmed temperature
  sim_yrs <- baseline_yrs + warming_yrs + final_yrs # total number of years
  Tmean <- baseT + c(rep(0,baseline_yrs),seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                     rep(deltaT,final_yrs))
  temperature = rnorm(sim_yrs,Tmean,Tstdev)  # generate temperature time series
  source("genetic_diversity_sim_temporal.R")
  
  # generate forecasts
  source("genetic_diversity_forecast.R")
  
  saved_weights[,iRun] <- weight
  
}

