
# call this script from genetic_diversity_master.R after running the previous scripts

forecast_start <- length(temporal_sim$temperature)-final_yrs-warming_yrs

forecast_Tmean <- Tmean[forecast_start:sim_yrs]

# start here spatial_forecast <- 