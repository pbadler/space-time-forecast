
# call this script from genetic_diversity_master.R

myCols<-c("blue","darkgrey","red","black")

# plot fecundity as function of temperature for both phenotypes
png("reactionnorms.png",height=3,width=4,res=300,units="in")
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  temp=seq(-8,8,0.01)
  react.norms=getLambdas(temp,fec.Tmu,fec.Tsigma,fec.max)
  matplot(temp,react.norms,type="l",xlab="Temperature",ylab="Germination rate",col=myCols,lty="solid")
  text(fec.Tmu[1],0.9*fec.max[1],"AA",col="blue")
  text(fec.Tmu[2],0.9*fec.max[2],"Aa",col="grey")
  text(fec.Tmu[3],0.9*fec.max[3],"aa",col="red")
dev.off()

# plot mean abundances vs mean temperature, along with spatial model "prediction"
xx <- seq(-5,5,0.01)
yy <- coef(spatial_model)[1] + coef(spatial_model)[2]*xx + coef(spatial_model)[3]*xx^2
png("spatial_pattern.png",height=3.5,width=5,units="in",res=400)
  par(mfrow=c(1,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
  matplot(Tmean,simN.mix[,c(2:5)],xlab="Mean temperature",ylab="N",
            type="l", col=myCols,
            lty=c("solid"),lwd=c(1,1,1,2))
  lines(xx,yy,lty="dashed",lwd=1)
  #legend("topright",c("AA","Aa","aa","Pop"),fill=myCols,bty="n")
dev.off()
