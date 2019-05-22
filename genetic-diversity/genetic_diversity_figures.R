
# call this script from genetic_diversity_master.R

myCols<-c("blue","darkgrey","red","black")

# figure to show spatial and temporal relationships
png("figures/spatial&temporal_model.png",height=10,width=4,res=400,units="in")
  
  par(mfrow=c(3,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1),cex=1)
  
  # plot fecundity as function of temperature for all genotypes
  temp=seq(-7,7,0.01)
  react_norms=getLambdas(temp,fec_Tmu,fec_Tsigma,fec_max)
  matplot(temp,react_norms,type="l",xlab="Temperature",ylab="Germination rate",col=myCols,lty="solid")
  text(fec_Tmu[1],0.9*fec_max[1],"AA",col=myCols[1])
  text(fec_Tmu[2],0.9*fec_max[2],"Aa",col=myCols[2])
  text(fec_Tmu[3],0.9*fec_max[3],"aa",col=myCols[3])
  mtext("(A)",side=3, line=0.5, adj=0)
  
  # plot mean abundances vs mean temperature, along with spatial model "prediction"
  xx <- seq(-5,5,0.01)
  yy <- coef(spatial_model)[1] + coef(spatial_model)[2]*xx + coef(spatial_model)[3]*xx^2
  matplot(simN_mix$Tmean,simN_mix[,c(2:5)],xlab="Mean temperature",ylab="N",
            type="l", col=myCols,
            lty=c("solid"),lwd=c(1,1,1,2))
  lines(xx,yy,lty="dashed",lwd=1)
  #legend("topright",c("AA","Aa","aa","Pop"),fill=myCols,bty="n")
  mtext("(B)",side=3, line=0.5, adj=0)
  
  # plot population growth and density vs annual temperature
  # first make colors based on N0 density
  rr <- range(N_t0_baseline)
  svals <- (N_t0_baseline-rr[1])/diff(rr)
  f <- colorRamp(c("red", "blue"), 0.4)
  density_colors <- alpha(rgb(f(svals)/255),0.6)
  # now plot
  plot(temperature_t1_baseline,r_baseline,xlab="Annual temperature",ylab="log Population growth rate",
       col=density_colors,pch=16)
  mtext("(C)",side=3, line=0.5, adj=0)

dev.off()


# plot forecast
png("figures/forecast.png",height=5,width=5,res=400,units="in")
  
  layout(matrix(c(1,2),2,1),heights=c(0.25,0.75),widths=1)
  
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(0,4,1,1))
  plot(temperature,xlab="",ylab="Temperature",type="l",xaxt="n",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines((baseline_yrs+1):sim_yrs,forecast_Tmean,lwd=2,col="black")
  
  par(mar=c(3,4,1,1))
  plot(N,xlab="Time",ylab="N",type="l",col="darkgrey")
  abline(v=baseline_yrs,col="black",lty=3)
  lines(1:length(N),spatial_forecast,col="red",lwd=2)
  lines(1:length(N),temporal_forecast,col="blue",lwd=2)
  lines(1:length(N),combined_forecast,col="purple",lwd=2,lty=2)
  legend("bottomleft",c("Observed","Spatial forecast","Temporal forecast","Combined forecast"),
         col=c("darkgrey","red","blue","purple"),lwd=c(1,2,2,2),lty=c(1,1,1,2),bty="n")

dev.off()

# plot weights and genotypes
png("figures/forecast_supplement.png",height=3.5,width=4.5,res=400,units="in")

  par(mfrow=c(1,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,4))
  
  matplot(N_genotypes[(baseline_yrs+1):sim_yrs,],type="l",xlab="Time",ylab="N",
          lty=1,col=myCols)
  legend(x=length(weight)*0.5,y=max(N_genotypes)*0.4,c("AA","Aa","aa","Weight"),
         lty=c(1,1,1,2),col=myCols,bty="n")
  
  par(new=T)
  plot(1:length(weight),weight,type="l",lty=2,xlab=NA,ylab=NA,axes=F)
  axis(side=4)
  mtext(side=4,line=3,"Weight")

dev.off()
