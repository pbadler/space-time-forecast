
# call this script from comm_turn_master.R

myCols<-c("blue3","red3")

# plot spatial and temporal models
png("figures/spatial_temporal_models.png",height=3.5,width=4,res=400,units="in")
  
plot(x=Tmean,y_spatial,type="n",xlab="Baseline temperature",ylab="Biomass")
points(x_temporal,y_temporal,pch=1,col=myCols[2])
abline(temporal_model,col=myCols[2], lwd=2)
points(x_spatial,y_spatial,pch=16,col=myCols[1])
abline(spatial_model,col=myCols[1],lwd=2)

dev.off()

# plot forecast
my_expected_temp <- Tmean[site] + c(rep(0,baseline_yrs),seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                          rep(deltaT,final_yrs))
my_obs_temp <- temperature[(burnin_yrs+1):sim_yrs,site]
my_biomass <- rowSums(spxp[(burnin_yrs+1):sim_yrs,site,])

png("figures/forecast.png",height=5,width=5,res=400,units="in")
  
  layout(matrix(c(1,2),2,1),heights=c(0.25,0.75),widths=1)
  
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(0,4,1,1))
  plot(my_obs_temp,xlab="",ylab="Temperature",type="l",xaxt="n",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(my_expected_temp,lwd=2,col="black")
  
  par(mar=c(3,4,1,1))
  plot(my_biomass,xlab="Time",ylab="N",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(1:length(my_biomass),spatial_forecast,col="red",lwd=2)
  lines(1:length(my_biomass),temporal_forecast,col="blue",lwd=2)
  lines(1:length(my_biomass),combined_forecast,col="purple",lwd=2,lty=2)
  legend("topleft",c("Observed","Spatial forecast","Temporal forecast","Combined forecast"),
         col=c("darkgrey","red","blue","purple"),lwd=c(1,2,2,2),lty=c(1,1,1,2),bty="n")

dev.off()

rm(my_biomass,my_obs_temp,my_expected_temp)

# mean temperature vs mean species biomass 
matplot(Tmean,spxp_mean_baseline,type="l",xlab="Baseline temperature",ylab="Biomass")

# biomass over time for all species at one location, baseline years
matplot(burnin_yrs:(burnin_yrs+baseline_yrs),
        spxp[burnin_yrs:(burnin_yrs+baseline_yrs),site,],
        type="l",xlab="Time",ylab="Biomass")

# total biomass at each location over time
matplot(1:sim_yrs,apply(spxp,MARGIN=c(1,2),FUN="sum"),type="l")

# # plot weights and genotypes
# png("figures/forecast_supplement.png",height=3.5,width=4.5,res=400,units="in")
# 
#   par(mfrow=c(1,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,4))
#   
#   matplot(N_genotypes[(baseline_yrs+1):sim_yrs,],type="l",xlab="Time",ylab="N",
#           lty=1,col=myCols)
#   legend(x=length(weight)*0.5,y=max(N_genotypes)*0.4,c("AA","Aa","aa","Weight"),
#          lty=c(1,1,1,2),col=myCols,bty="n")
#   
#   par(new=T)
#   plot(1:length(weight),weight,type="l",lty=2,xlab=NA,ylab=NA,axes=F)
#   axis(side=4)
#   mtext(side=4,line=3,"Weight")
# 
# dev.off()
