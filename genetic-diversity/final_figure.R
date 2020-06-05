
# call this script from genetic_diversity_master.R

###
### 1. Fig. 6 
###

# read in weights from community turnover example
saved_weights_comm <- read.csv("./../community-turnover/simulations/saved_weights.csv",header=T)
d_vals <- colnames(saved_weights_comm)
d_vals <- gsub("X","",d_vals)

# Figure to compare rate of change in weights for different
# widths of temperature niche
myCols<-c("black","gray40","gray80")
png("figures/dispersal_niche_width.png",height=3.5,width=8.5,res=400,units="in")
  
  par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))
  
  #community example
  matplot(saved_weights_comm,type="l",lwd=2,lty=1,col=myCols,
          xlab="Time step",ylab="Weight")
  legend("topright",legend=d_vals,
       title="dispersal",lwd=2,lty=1,col=myCols,bty="n")
  mtext("(A)",side=3,line=0.5,adj=0)
  
  #eco-evo example
  matplot(saved_weights,type="l",lwd=2,lty=1,col=myCols,
          xlab="Time step",ylab="Weight")
  legend("topright",legend=fec_Tsigma_baseline[1]*delta_fec_sigma,
         title=expression(sigma[T]),lwd=2,lty=1,col=myCols,bty="n")
  mtext("(B)",side=3,line=0.5,adj=0)
  
dev.off()


###
### 2. Supp. Fig on variation among replicate runs
###

# read in weights from community turnover example
saved_weights_comm <- read.csv("./../community-turnover/simulations/saved_weights_reps.csv",header=T)

png("figures/replicates_weights.png",height=3.5,width=8.5,res=400,units="in")

par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

#community example
matplot(saved_weights_comm,type="l",lwd=2,lty=1,col="black",
        xlab="Time step",ylab="Weight")
mtext("(A)",side=3,line=0.5,adj=0)

#eco-evo example
matplot(saved_weights_reps,type="l",lwd=2,lty=1,col="black",
        xlab="Time step",ylab="Weight")
mtext("(B)",side=3,line=0.5,adj=0)

dev.off()
