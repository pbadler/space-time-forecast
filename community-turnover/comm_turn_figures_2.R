

# call this script from comm_turn_master.R

# Figure to compare rate of change in weights for different
# dispersal rates
myCols<-c("black","gray40","gray80")
png("figures/compare_dispersal.png",height=3.5,width=5,res=400,units="in")
par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1))
matplot(saved_weights_d,type="l",lwd=2,lty=1,col=myCols,
        xlab="Time step",ylab="Weight")
legend("topright",legend=d_vals,
       title="dispersal",lwd=2,lty=1,col=myCols,bty="n")
dev.off()
