# call this script from comm_turn_master.R

# set up landscape
landscape <- cbind(1:L_land , rep(1 ,L_land),Tmean)

# 1st step: initialisation of the metacommunity
spxp <- array(NA, dim=c(L_land, N, sim_yrs))
spxp[,,1] <- 2/N
dt <- 0.025   # ????
d <- d/dt

# Generation of the species pool
tr <- SpeciesPoolGen(N, Tmin, Tmax, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)

# dispersal matrix (seeds disperse only to nearest neighbor locations)
seed_rain <- as.matrix(dist(landscape[,1:2],diag=T,upper=T))
seed_rain[seed_rain > 1] <- 0
seed_rain <- seed_rain*(d*dt*0.5)
diag(seed_rain) <- 1-d*dt
seed_rain[1,1] = 1-dt*d/2  
seed_rain[L_land,L_land] = 1-dt*d/2

# Run of the initialisation of community
for(t in 2:sim_yrs){
  yr_temp <- temperature[t,]
  spxp[,,t] <- CommunityTempDis(spxp[,,t-1], tr, seed_rain, yr_temp, dt)
} 

spxp_baseline <- apply(spxp[,,1:baseline_yrs],MARGIN = c(1,2),FUN=mean)
# save the initial metacommunity
dir.create("simulations/run_1")
write.csv(spxp,"simulations/run_1/spxp.csv",row.names = F)
write.csv(landscape,"simulations/run_1/landscape.csv",row.names = F)
write.csv(tr,"simulations/run_1/tr.csv",row.names = F)

