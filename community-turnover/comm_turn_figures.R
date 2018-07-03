
# call this script from comm_turn_master.R

myCols<-alpha(c("blue3","red3"),0.4)

###
### Total biomass figures
###

# plot spatial and temporal models
png("figures/community_models.png",height=3.5,width=4,res=400,units="in")

  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  
  # set up empty axes
  plot(x=Tmean,y_spatial,type="n",xlab="Temperature",ylab="Biomass")
  abline(v=Tmean[site], lty=2,col="darkgray")  # show focal site
  
  # total biomass temporal
  xx <- seq(min(x_temporal),max(x_temporal),0.1)
  points(x_temporal,y_temporal,pch=16,col=myCols[2])
  y_hat <- coef(temporal_model)[1] + coef(temporal_model)[2]*mean(y_temporal) + coef(temporal_model)[3]*xx
  lines(xx,y_hat, col=myCols[2], lwd=2 ) # plot marginal effect of temperature for temporal model
  
  # focal species temporal
  points(x_temporal,y_temporal_spp,pch=1,col=myCols[2])
  y_hat <- coef(temporal_model_spp)[1] + coef(temporal_model_spp)[2]*mean(y_temporal_spp) + coef(temporal_model_spp)[3]*xx
  lines(xx,y_hat, col=myCols[2], lwd=1, lty="dotted" ) # plot marginal effect of temperature for temporal model
  
  # total biomass spatial
  points(x_spatial,y_spatial,pch=16,col=myCols[1])
  abline(spatial_model,col=myCols[1],lwd=2)
  
  # focal species spatial
  points(x_spatial,y_spatial_spp,pch=1,col=myCols[1])
  y_hat <-  coef(spatial_model_spp)[1] + coef(spatial_model_spp)[2]*x_spatial + coef(spatial_model_spp)[3]*x_spatial^2
  lines(x_spatial,y_hat,col=myCols[1],lwd=1, lty="dotted")
  
  legend("topleft",c("Spatial (total)","Temporal (total)","Spatial (species)","Temporal (species)"),
         lty=c("solid","solid","dotted","dotted"),pch=c(16,16,1,1),lwd=c(2,2,1,1),
         col=c(myCols,myCols),bty="n",cex=0.8)

dev.off()


# plot forecast for total biomass
my_expected_temp <- Tmean[site] + c(rep(0,baseline_yrs),seq(deltaT/warming_yrs,deltaT,deltaT/warming_yrs),
                          rep(deltaT,final_yrs))
my_obs_temp <- temperature[(burnin_yrs+1):sim_yrs,site]
my_biomass <- rowSums(spxp[(burnin_yrs+1):sim_yrs,site,])
my_biomass_spp <- spxp[(burnin_yrs+1):sim_yrs,site,my_species]

png("figures/community_forecast.png",height=7.5,width=4,res=400,units="in")
  
  layout(matrix(c(1,2,3),3,1),heights=c(0.1,0.4,0.4),widths=1)
  
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(1,4,1,1),oma=c(2,0,0,0),cex=0.8)
  
  plot(my_obs_temp,xlab="",ylab="Temperature",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(my_expected_temp,lwd=2,col="black")
  
  # total biomass
  plot(my_biomass,xlab="",ylab="Total biomass",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(1:length(my_biomass),spatial_forecast,col="red",lwd=2)
  lines(1:length(my_biomass),temporal_forecast,col="blue",lwd=2)
  lines(1:length(my_biomass),combined_forecast,col="purple",lwd=2,lty=2)

  # focal species
  plot(my_biomass_spp,xlab="",ylab="Focal species biomass",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(1:length(my_biomass_spp),spatial_forecast_spp,col="red",lwd=2)
  lines(1:length(my_biomass_spp),temporal_forecast_spp,col="blue",lwd=2)
  lines(1:length(my_biomass_spp),combined_forecast_spp,col="purple",lwd=2,lty=2)
  
  legend("right",c("Observed","Spatial forecast","Temporal forecast","Combined forecast"),
         col=c("darkgrey","red","blue","purple"),lwd=c(1,2,2,2),lty=c(1,1,1,2),bty="n",cex=0.8)
  
  mtext(side=1,"Time",outer=T,line=1,cex=0.8)
  
dev.off()

rm(my_biomass,my_biomass_spp,my_obs_temp,my_expected_temp)


# plot community change and weights over time 

png("figures/community_change_plus_weights.png",height=3.5,width=4.5,units="in",res=400)

  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,4))

  spp_colors <- rep("darkgrey",N)
  spp_colors[my_species] <- "black"
  matplot(burnin_yrs:sim_yrs,
          spxp[burnin_yrs:sim_yrs,site,],
          type="l",lty=1,xlab="Time",ylab="Biomass by species",col=spp_colors)
  abline(v=burnin_yrs+baseline_yrs+1,lty="dashed")
  
  # add weights
  par(new=T)
  plot((burnin_yrs+1):sim_yrs,c(rep(NA,baseline_yrs),weight),type="l",lty=1,lwd=2,col="dodgerblue",
       xlab=NA,ylab=NA,ylim=c(0,1),axes=F)
  axis(side=4)
  mtext(side=4,line=2,"Weight")
  lines((burnin_yrs+1):sim_yrs,c(rep(NA,baseline_yrs),weight_spp),lty=1,lwd=2,col="firebrick1")

  legend("topright",c("Total","Species"),lwd=2,col=c("dodgerblue","firebrick1"),
         bty="n",title="Weights",cex=0.8)
  
dev.off()


# mean temperature vs mean species biomass 
png("figures/mean_biomass_spp_by_site.png",height=3.4,width=4,res=400,units="in")
  
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,4))
  y <-spxp_mean_baseline[,colSums(spxp_mean_baseline)>0]
  jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"),bias=2.5)
  tmp <- jet.colors(n=NCOL(y))
  matplot(Tmean,y,type="l",xlab="Baseline temperature",ylab="Biomass",
          lty=1,col=tmp)

dev.off()

# total biomass at each location over time
matplot(1:sim_yrs,apply(spxp,MARGIN=c(1,2),FUN="sum"),type="l")

# # plot weights vs community change 

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

