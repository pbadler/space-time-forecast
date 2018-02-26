
# clean up
rm(list=ls())

# set working directory
setwd("C:/Repos/space-time-forecast/genetic-diversity")

# load packages and functions
source("genetic_diversity_functions.R")


###
### 1. Describe spatial pattern in abundances across a range of mean temperatures
###

sim.yrs<-2000
Tmean = seq(-5,5,0.25)   # set of mean temperatures to explore
Tstdev <- 1   # st dev of temperature

# generate one set of temperatures to use in all simulations
temperature<-matrix(NA,sim.yrs,length(Tmean))
for(iT in 1:length(Tmean)){
  temperature[,iT] = rnorm(sim.yrs,Tmean[iT],Tstdev)  # generate temperature time series
}  

# diploid model
simName <- "ss.5_noDominant"
fec.Tmu = c(-1,0,1)  # optimal temperature for genotypes AA, Aa, and aa
fec.Tsigma = rep(5,3)     # standard deviation in temperature response for genotypes AA, Aa, and aa
fec.max = c(100,100,100)  # fecundity for genotypes AA, Aa, and aa
G <- 1
seedSurv = 0.5  # survival of ungerminated seeds (same for both genotypes)
alpha = matrix(1,3,3)      # all competition coefficients (intra- and inter-) are equal for both genotypes

source("genetic_diversity_experiment.R")

#CHECK: why zeros for least fit genotype's mean abundances at lowest temps but not highest?

# FIGURES

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
