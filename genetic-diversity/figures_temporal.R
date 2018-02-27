
# show key features of temporal model
png("figures/temporal_pattern.png",height=4,width=7.5,res=400,units="in")
  par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  # plot temperature vs population growth during baseline period
  plot(temperature_t1_baseline,r_baseline,xlab="Temperature",ylab="log Population growth rate")
  # plot observed vs predicted population growth
  plot(r_baseline,predict(temporal_model),xlab="Observed", ylab= "Predicted")
dev.off()


# plot time series with spatial and temporal predictions!
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