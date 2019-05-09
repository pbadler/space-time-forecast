
# call this script from genetic_diversity_master.R

# Figure to compare rate of change in weights for different
# widths of temperature niche
myCols<-c("black","gray40","gray80")
png("figures/compare_Tsigmas.png",height=3.5,width=5,res=400,units="in")
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  matplot(saved_weights,type="l",lwd=2,lty=1,col=myCols,
          xlab="Time step",ylab="Weight")
  legend("topright",legend=fec_Tsigma_baseline[1]*delta_fec_sigma,
         title=expression(sigma[T]),lwd=2,lty=1,col=myCols,bty="n")
dev.off()


