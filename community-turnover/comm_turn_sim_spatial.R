# call this script from comm_turn_master.R

# 1st step: initialisation of the metacommunity
dt <- 0.025
spxp <- matrix(1,L_land, N)*(2/N)
d <- d/dt

# Generation of the species pool
tr <- SpeciesPoolGen(N, Tmin, Tmax, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)
temp <- matrix(landscape[,3],L_land,1)

# dispersal matrix (seeds disperse only to nearest neighbor locations)
seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
seed_rain[seed_rain > 1] <- 0
seed_rain <- seed_rain*(d*dt*0.5)
diag(seed_rain) <- 1-d*dt
seed_rain[1,1] = 1-dt*d/2  
seed_rain[L_land,L_land] = 1-dt*d/2

# Run of the initialisation of community
for(t in 1:sim_yrs){
  spxp <- CommunityTempDis(spxp, tr, seed_rain, temp, dt)
} 
spxp_ini <- spxp

# save the initial metacommunity
dir.create("simulations/run_1")
write.csv(spxp_ini,"simulations/run_1/spxp_ini.csv",row.names = F)
write.csv(landscape,"simulations/run_1/landscape.csv",row.names = F)
write.csv(tr,"simulations/run_1/tr.csv",row.names = F)

