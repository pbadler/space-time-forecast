
# call this script from genetic_diversity_master.R

myCols<-c("blue","darkgrey","red","black")

# plot fecundity as function of temperature for both phenotypes
png("figures/reactionnorms.png",height=3,width=4,res=300,units="in")
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  temp=seq(-8,8,0.01)
  react_norms=getLambdas(temp,fec_Tmu,fec_Tsigma,fec_max)
  matplot(temp,react_norms,type="l",xlab="Temperature",ylab="Germination rate",col=myCols,lty="solid")
  text(fec_Tmu[1],0.9*fec_max[1],"AA",col="blue")
  text(fec_Tmu[2],0.9*fec_max[2],"Aa",col="grey")
  text(fec_Tmu[3],0.9*fec_max[3],"aa",col="red")
dev.off()

# plot mean abundances vs mean temperature, along with spatial model "prediction"
xx <- seq(-5,5,0.01)
yy <- coef(spatial_model)[1] + coef(spatial_model)[2]*xx + coef(spatial_model)[3]*xx^2
png("figures/spatial_pattern.png",height=3.5,width=5,units="in",res=400)
  par(mfrow=c(1,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  matplot(Tmean,simN_mix[,c(2:5)],xlab="Mean temperature",ylab="N",
            type="l", col=myCols,
            lty=c("solid"),lwd=c(1,1,1,2))
  lines(xx,yy,lty="dashed",lwd=1)
  #legend("topright",c("AA","Aa","aa","Pop"),fill=myCols,bty="n")
dev.off()
