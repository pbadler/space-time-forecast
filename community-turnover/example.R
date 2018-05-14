#
#
#	Example script to run the model of Alexander et al. 2018 'Lags in the response of mountain plant communities to climate change' - Global Change Biology
#	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu
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
tr = SpeciesPoolGen(N, 4, 16.5, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)
temp = matrix(landscape[,3],L_land,1)

# dispersal matrix
# original matlab code: seed_rain =  squareform(pdist(landscape(:,1:2)) < 2).*(d*dt*0.5);
seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))*(d*dt*0.5)  # where does the < 2 go?
indx <-  (diag(rep(1,dim(landscape)[1])) == 1)
seed_rain[indx==T] <- 1-d*dt
seed_rain[1,1] = 1-dt*d/2  
seed_rain[dim(landscape)[1],dim(landscape)[1]] = 1-dt*d/2

# Run of the initialisation of community
for(t in 1:Tmax){
  spxp <- CommunityTempDis(spxp, tr, seed_rain, temp, dt)
} 
spxp_ini <- spxp

# save the initial metacommunity
dir.create("simulations/run_1")
write.csv(spxp_ini,"simulations/run_1/spxp_ini.txt",row.names = F)
write.csv(landscape,"simulations/run_1/landscape.txt",row.names = F)
write.csv(tr,"simulations/run_1/tr.txt",row.names = F)

spxp_ini[spxp_ini < 10e-12] = 0 # apply extinction threshold

spxp <- spxp_ini
# From the initial situation, re-run community dynamic while temperature is increasing.
temp <- landscape[,3]
time_to_max <- 20000
temp_increase <- 3
for(t in 1:time_to_max){
  temp_new <- temp_increase*t/time_to_max + temp
  spxp <- CommunityTempDis(spxp, tr, seed_rain, temp_new, dt)
}
spxp_final <- spxp



