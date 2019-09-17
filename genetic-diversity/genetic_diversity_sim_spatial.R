
# simulate heterozygous population
# most parameters are assigned in genetic_diversity_master.R

burn_in=sim_yrs/5   # years to reach stochastic equilibrium

simN_mix <- data.frame("Tmean"=0,"N_AA"=0,"N_Aa"=0,"N_aa"=0,"N_Pop"=0)
counter=0
for(iT in 1:length(Tmean)){
  
  # initialize
  N=Plants=Fec=matrix(0,sim_yrs,3)
  colnames(N)<-colnames(Plants)<-colnames(Fec)<-c("AA","Aa","aa")
  N[1,] = fec_max[1]/3
  lambdas = getLambdas(temperature[,iT],fec_Tmu,fec_Tsigma,fec_max)
  
  for(iYr in 2:(sim_yrs)){
    seeds<-newSeeds(N[iYr-1,],lambdas[iYr,],alpha,G)
    phi <- getPhi(seeds)
    out <- diploid(N[iYr-1,],seedSurv,G,seeds,phi)
    N[iYr,] <- out$SB
    Plants[iYr,] <- out$Plants
    Fec[iYr,] <- out$SP
  }   # next i
  
  myOut<-processOutput(N,burn_in,temperature[,iT])
  # save results
  
  counter <- counter + 1
  simN_mix[counter,] <-c(Tmean[iT],myOut$means)

}

simN_mix[is.na(simN_mix)] <- 0 # replace missing values with zero

# fit a "spatial" model describing mean abundance as a function of mean temperature
Tmean2 <- simN_mix$Tmean^2
spatial_model <- lm(N_Pop~Tmean + Tmean2,data=simN_mix)
