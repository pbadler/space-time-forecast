
# call this script from comm_turn_master.R after running the previous scripts

###
### Model total biomass
###

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
# make predictions from spatial model with uncertainty
# spatial_forecast_spp <- predict(spatial_model_spp,newdata=data.frame(x_spatial=forecast_Tmean),
#                                 interval="predict")

# make predictions from temporal model
temporal_forecast <- rep(NA,length(forecast_Tmean)+1)
temporal_forecast[1] <- y_temporal[length(y_temporal)] # initialize at biomass from last baseline year
for(iTime in 1:length(forecast_Tmean)){
   temporal_forecast[iTime+1] <- coef(temporal_model)[1] +coef(temporal_model)[2]*temporal_forecast[iTime] + 
                                  coef(temporal_model)[3]*forecast_Tmean[iTime]
}
temporal_forecast <- temporal_forecast[2:length(temporal_forecast)]  # drop first ("observed") year

# make predictions from temporal model using a Monte Carlo approach for uncertainty
# ndraws <- 100
# pars <- rmvnorm(ndraws,coef(temporal_model_spp),vcov(temporal_model_spp))
# temporal_forecast_spp <- matrix(NA,length(forecast_Tmean)+1,ndraws)
# temporal_forecast_spp[1,] <- y_temporal_spp[length(y_temporal)] # initialize at biomass from last baseline year
# for(iPar in 1:ndraws){
#   for(iTime in 1:length(forecast_Tmean)){
#     temporal_forecast_spp[iTime+1,iPar] <- pars[iPar,1] +pars[iPar,2]*temporal_forecast_spp[iTime,iPar] + 
#       pars[iPar,3]*forecast_Tmean[iTime]
#   }
# }
# temporal_forecast_spp <- temporal_forecast_spp[2:dim(temporal_forecast_spp)[1],]  # drop first ("observed") year

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

###
###  Model biomass of the dominant species 
###

my_species <- which(spxp_mean_baseline[site,]==max(spxp_mean_baseline[site,]))

# spatial model
x_spatial <- Tmean
y_spatial_spp <- spxp_mean_baseline[,my_species]
spatial_model_spp <- lm(y_spatial_spp~x_spatial+I(x_spatial^2))

# temporal model
site <- L_land/2   # focal site
x_temporal <- temperature[burnin_yrs:(burnin_yrs+baseline_yrs),site]
y_temporal_spp <- spxp[burnin_yrs:(burnin_yrs+baseline_yrs),site,my_species]
y_lag_spp <- spxp[(burnin_yrs-1):(burnin_yrs+baseline_yrs-1),site,my_species]
temporal_model_spp <- lm(y_temporal_spp ~ y_lag_spp + x_temporal)

# get mean temperatures for forecast period (after baseline_yrs)
forecast_Tmean <- Tmean[site] + c(seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                                  rep(deltaT,final_yrs))

# make predictions from spatial model
spatial_forecast_spp <- coef(spatial_model_spp)[1] + coef(spatial_model_spp)[2]*forecast_Tmean +
                            coef(spatial_model_spp)[3]*forecast_Tmean^2

# make predictions from temporal model
temporal_forecast_spp <- rep(NA,length(forecast_Tmean)+1)
temporal_forecast_spp[1] <- y_temporal_spp[length(y_temporal)] # initialize at biomass from last baseline year
for(iTime in 1:length(forecast_Tmean)){
  temporal_forecast_spp[iTime+1] <- coef(temporal_model_spp)[1] +coef(temporal_model_spp)[2]*temporal_forecast_spp[iTime] + 
    coef(temporal_model_spp)[3]*forecast_Tmean[iTime]
}
temporal_forecast_spp <- temporal_forecast_spp[2:length(temporal_forecast_spp)]  # drop first ("observed") year

# fit combined model
y <- spxp[(burnin_yrs+baseline_yrs+1):sim_yrs,site,my_species]
betas<-c(10,-1) # initial values for regression coefficients
forecast_yrs <- warming_yrs+final_yrs
predict_error <- function(betas){
  weight <- inv.logit(betas[1] + betas[2]*(sqrt(1:forecast_yrs)))
  y_hat <- weight*temporal_forecast_spp + (1-weight)*spatial_forecast_spp
  sum_sq_errors <- sum((y-y_hat)^2)
  return(sum_sq_errors)
}
fitted <- optim(betas,predict_error)

# make forecast with combined model
weight_spp <- inv.logit(fitted$par[1] + fitted$par[2]*(sqrt(1:forecast_yrs)))
combined_forecast_spp <- weight_spp*temporal_forecast_spp + (1-weight_spp)*spatial_forecast_spp

# add NAs to forecast time series so they line up with simulated ("observed") time series
spatial_forecast_spp <- c(rep(NA,baseline_yrs),spatial_forecast_spp)
temporal_forecast_spp <- c(rep(NA,baseline_yrs),temporal_forecast_spp)
combined_forecast_spp <-c(rep(NA,baseline_yrs),combined_forecast_spp)
