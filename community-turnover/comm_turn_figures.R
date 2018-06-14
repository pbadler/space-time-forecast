
# call this script from comm_turn_master.R

myCols<-c("blue3","red3")

# plot spatial and temporal models
png("figures/spatial_temporal_models.png",height=3.5,width=4,res=400,units="in")

  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  
  plot(x=Tmean,y_spatial,type="n",xlab="Temperature",ylab="Biomass")
  points(x_temporal,y_temporal,pch=1,col=myCols[2])
  y_hat <- coef(temporal_model)[1] + coef(temporal_model)[2]*mean(y_temporal) + coef(temporal_model)[3]*x_temporal
  lines(x_temporal,y_hat, col=myCols[2], lwd=2 ) # plot marginal effect of temperature for temporal model
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
  plot(my_biomass,xlab="Time",ylab="Biomass",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(1:length(my_biomass),spatial_forecast,col="red",lwd=2)
  lines(1:length(my_biomass),temporal_forecast,col="blue",lwd=2)
  lines(1:length(my_biomass),combined_forecast,col="purple",lwd=2,lty=2)
  legend("topleft",c("Observed","Spatial forecast","Temporal forecast","Combined forecast"),
         col=c("darkgrey","red","blue","purple"),lwd=c(1,2,2,2),lty=c(1,1,1,2),bty="n")

dev.off()

rm(my_biomass,my_obs_temp,my_expected_temp)

# plot community change and weights over time 

png("figures/community_change_plus_weight.png",height=3.5,width=4.5,units="in",res=400)

  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,4))

  matplot(burnin_yrs:sim_yrs,
          spxp[burnin_yrs:sim_yrs,site,],
          type="l",xlab="Time",ylab="Biomass")
  abline(v=burnin_yrs+baseline_yrs+1,lty="dashed")
  
  # add weights
  par(new=T)
  plot((burnin_yrs+1):sim_yrs,c(rep(NA,baseline_yrs),weight),type="l",lty=1,lwd=2,col="darkgray",
       xlab=NA,ylab=NA,axes=F)
  axis(side=4)
  mtext(side=4,line=2,"Weight")

dev.off()

# mean temperature vs mean species biomass 
matplot(Tmean,spxp_mean_baseline,type="l",xlab="Baseline temperature",ylab="Biomass")

# total biomass at each location over time
matplot(1:sim_yrs,apply(spxp,MARGIN=c(1,2),FUN="sum"),type="l")



# # plot weights vs turnover

# calculate Euclidean distance 
comm_matrix <- spxp[(burnin_yrs+baseline_yrs):sim_yrs,site,]

# Euclidean distance
ED = function(x,y){
  #x and y are vectors of spp abundances
  #they must be the same length!
  if(length(x)!=length(y)) stop("Bad abundances!")
  out =  sqrt(sum((x-y)^2))
  out
}

obs_ED <- rep(NA,warming_yrs+final_yrs)
for(i in 1:length(obs_ED)){
  obs_ED[i] <- ED(x=comm_matrix[1,],y=comm_matrix[1+i,])
}

jaccards=function(x,y){
  #x and y are vectors of spp abundances
  #they must be the same length!
  if(length(x)!=length(y)) stop("Bad abundances!")
  # get number of species present at both sites
  a = length(which(x>0 & y>0))
  # get number of species present at first but not second site
  b = length(which(x>0 & y==0))
  # get number of species present at second but not first site
  c = length(which(x==0 & y>0))
  out = a/(a+b+c)
  out
}

obs_JC <- rep(NA,warming_yrs+final_yrs)
for(i in 1:length(obs_JC)){
  obs_JC[i] <- jaccards(x=comm_matrix[1,],y=comm_matrix[1+i,])
}


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
