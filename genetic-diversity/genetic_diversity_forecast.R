
# call this script from genetic_diversity_master.R after running the previous scripts

# get mean temperatures for forecast period (after baseline_yrs)
tmp <- length(Tmean)-final_yrs-warming_yrs
forecast_Tmean <-  Tmean[(tmp+1):length(Tmean)]

# make predictions from spatial model
spatial_forecast <- coef(spatial_model)[1] + coef(spatial_model)[2]*forecast_Tmean + coef(spatial_model)[3]*forecast_Tmean^2

# make predictions from temporal model
temporal_forecast <- rep(NA,length(forecast_Tmean)+1)
temporal_forecast[1] <- N[baseline_yrs] # initialize at N from last baseline year
for(iTime in 1:length(forecast_Tmean)){
  r_iTime <- coef(temporal_model)[1] +coef(temporal_model)[2]*log(temporal_forecast[iTime]) + 
                                  coef(temporal_model)[3]*forecast_Tmean[iTime] + coef(temporal_model)[4]*forecast_Tmean[iTime]^2
  temporal_forecast[iTime+1] <- exp(r_iTime)*temporal_forecast[iTime]
}
temporal_forecast <- temporal_forecast[2:length(temporal_forecast)]  # drop first ("observed") year

# fit combined model
y <- N[(baseline_yrs+1):sim_yrs]
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


