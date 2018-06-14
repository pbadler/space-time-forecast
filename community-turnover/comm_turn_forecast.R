
# call this script from comm_turn_master.R after running the previous scripts

# spatial model
x_spatial <- Tmean
y_spatial <- rowSums(spxp_mean_baseline)
spatial_model <- lm(y_spatial~x_spatial)

# temporal model
site <- L_land/2   # focal site
x_temporal <- temperature[burnin_yrs:(burnin_yrs+baseline_yrs),site]
y_temporal <- rowSums(spxp[burnin_yrs:(burnin_yrs+baseline_yrs),site,])
y_lag <- rowSums(spxp[(burnin_yrs-1):(burnin_yrs+baseline_yrs-1),site,])
temporal_model <- lm(y_temporal ~ y_lag + x_temporal)

# get mean temperatures for forecast period (after baseline_yrs)
forecast_Tmean <- Tmean[site] + c(seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                          rep(deltaT,final_yrs))

# make predictions from spatial model
spatial_forecast <- coef(spatial_model)[1] + coef(spatial_model)[2]*forecast_Tmean 

# make predictions from temporal model
temporal_forecast <- rep(NA,length(forecast_Tmean)+1)
temporal_forecast[1] <- y_temporal[length(y_temporal)] # initialize at biomass from last baseline year
for(iTime in 1:length(forecast_Tmean)){
   temporal_forecast[iTime+1] <- coef(temporal_model)[1] +coef(temporal_model)[2]*temporal_forecast[iTime] + 
                                  coef(temporal_model)[3]*forecast_Tmean[iTime]
}
temporal_forecast <- temporal_forecast[2:length(temporal_forecast)]  # drop first ("observed") year


# fit combined model
y <- rowSums(spxp[(burnin_yrs+baseline_yrs+1):sim_yrs,site,])
betas<-c(10,-1) # initial values for regression coefficients
forecast_yrs <- warming_yrs+final_yrs
predict_error <- function(betas){
  weight <- inv.logit(betas[1] + betas[2]*(sqrt(1:forecast_yrs)))
  y_hat <- weight*temporal_forecast + (1-weight)*spatial_forecast
  sum_sq_errors <- sum((y-y_hat)^2)
  return(sum_sq_errors)
}
fitted <- optim(betas,predict_error)

# make forecast with combined model
weight <- inv.logit(fitted$par[1] + fitted$par[2]*(sqrt(1:forecast_yrs)))
combined_forecast <- weight*temporal_forecast + (1-weight)*spatial_forecast

# add NAs to forecast time series so they line up with simulated ("observed") time series
spatial_forecast <- c(rep(NA,baseline_yrs),spatial_forecast)
temporal_forecast <- c(rep(NA,baseline_yrs),temporal_forecast)
combined_forecast <-c(rep(NA,baseline_yrs),combined_forecast)


