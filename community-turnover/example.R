#
#
#	Example script to run the model of Alexander et al. 2018 'Lags in the response of mountain plant communities to climate change' - Global Change Biology
#	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu, https://github.com/LoicChr
# Translated to R by Peter Adler - peter.adler@usu.edu

rm(list=ls())

# set working directory
setwd("C:/Repos/space-time-forecast/community-turnover")

# read in functions
source("lib/CommunityTempDis.R")
source("lib/SpeciesPoolGen.R")

# generation of the Landscape aka a table with the following characteristics: 
# X coordinate, Y coordinate, temperature
L_land <- 10
landscape <- cbind(1:L_land , rep(1 ,L_land), seq(0,20, length=L_land))

# Parameters to generate the species pool
N <- 5
Gmax <- 0.5
Gmin <- 0.2
Lmax <- 1.5
Lmin <- 0.7
Cmax <- 0.2
Cmin <- 0.2 
d <- 10^(-3)

# 1st step: initialisation of the metacommunity
Tmax <- 40000; dt <- 0.025
spxp <- matrix(1,L_land, N)*(2/N)
d <- d/dt

# Generation of the species pool
tr <- SpeciesPoolGen(N, 4, 16.5, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)
temp <- matrix(landscape[,3],L_land,1)

# dispersal matrix (seeds disperse only to nearest neighbor locations)
seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
seed_rain[seed_rain > 1] <- 0
seed_rain <- seed_rain*(d*dt*0.5)
diag(seed_rain) <- 1-d*dt
seed_rain[1,1] = 1-dt*d/2  
seed_rain[L_land,L_land] = 1-dt*d/2

# Run of the initialisation of community
for(t in 1:Tmax){
  spxp <- CommunityTempDis(spxp, tr, seed_rain, temp, dt)
} 
spxp_ini <- spxp

# save the initial metacommunity
dir.create("simulations/run_1")
write.csv(spxp_ini,"simulations/run_1/spxp_ini.csv",row.names = F)
write.csv(landscape,"simulations/run_1/landscape.csv",row.names = F)
write.csv(tr,"simulations/run_1/tr.csv",row.names = F)

spxp_ini[spxp_ini < 10e-12] = 0 # apply extinction threshold

spxp <- spxp_ini
# From the initial situation, re-run community dynamic while temperature is increasing.
#temp <-  matrix(landscape[,3],L_land,1) # redundant?
time_to_max <- 20000
temp_increase <- 3
for(t in 1:time_to_max){
  temp_new <- temp_increase*t/time_to_max + temp
  spxp <- CommunityTempDis(spxp, tr, seed_rain, temp_new, dt)
}
spxp_final <- spxp

# figures
par(mfrow=c(1,2),tcl=-0.1,mgp=c(2,0.5,0),mar=c(3,4,3,1))
max_bio <- max(max(spxp_ini),spxp_final)
matplot(x=landscape[,3],y=spxp_ini,type="l",lty=1,ylim=c(0,max_bio*1.1),xlab="Temperature",ylab="Biomass")
mtext(side=3,"Initial",adj=0)
matplot(x=landscape[,3]+temp_increase,y=spxp_final,type="l",lty=1,ylim=c(0,max_bio*1.02),xlab="Temperature",ylab="Biomass")
mtext(side=3,"Final",adj=0)

# total biomass
plot(x=landscape[,3],rowSums(spxp_ini))
