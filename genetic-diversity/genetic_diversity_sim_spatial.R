
# simulate heterozygous population
# most parameters are assigned in genetic_diversity_master.R

burn.in=sim.yrs/5   # years to reach stochastic equilibrium

simN.mix <- data.frame("Tmean"=0,"N_AA"=0,"N_Aa"=0,"N_aa"=0,"N_Pop"=0)
counter=0
for(iT in 1:length(Tmean)){
  
  # initialize
  N=Plants=Fec=matrix(0,sim.yrs,3)
  colnames(N)<-colnames(Plants)<-colnames(Fec)<-c("AA","Aa","aa")
  N[1,] = fec.max[1]/3
  lambdas = getLambdas(temperature[,iT],fec.Tmu,fec.Tsigma,fec.max)
  
  for(iYr in 2:(sim.yrs)){
    seeds<-newSeeds(N[iYr-1,],lambdas[iYr,],alpha,G)
    phi <- getPhi(seeds)
    out <- diploid(N[iYr-1,],seedSurv,G,seeds,phi)
    N[iYr,] <- out$SB
    Plants[iYr,] <- out$Plants
    Fec[iYr,] <- out$SP
  }   # next i
  
  myOut<-processOutput(N,burn.in,temperature[,iT])
  # save results
  
  counter <- counter + 1
  simN.mix[counter,] <-c(Tmean[iT],myOut$means)

}

simN.mix[is.na(simN.mix)] <- 0 # replace missing values with zero

# fit a "spatial" model describing mean abundance as a function of mean temperature
spatial_model <- lm(N_Pop~Tmean + I(Tmean^2),data=simN.mix)
