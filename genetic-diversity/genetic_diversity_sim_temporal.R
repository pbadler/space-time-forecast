
# simulate a heterozygous population
# parameters are assigned in genetic_diversity_master.R

burn_in=baseline_yrs/2   # years to reach stochastic equilibrium
  
# initialize
N=Plants=Fec=matrix(0,sim_yrs,3)
colnames(N)<-colnames(Plants)<-colnames(Fec)<-c("AA","Aa","aa")
N[1,] = fec.max[1]/3
lambdas = getLambdas(temperature,fec.Tmu,fec.Tsigma,fec.max)

for(iYr in 2:(sim_yrs)){
  seeds<-newSeeds(N[iYr-1,],lambdas[iYr,],alpha,G)
  phi <- getPhi(seeds)
  out <- diploid(N[iYr-1,],seedSurv,G,seeds,phi)
  N[iYr,] <- out$SB
  Plants[iYr,] <- out$Plants
  Fec[iYr,] <- out$SP
}   # next i

# calculate per capita growth rates
temporal_sim<-processOutput(N,burn_in,temperature)

# fit a "temporal" model describing mean abundance as a function of mean temperature
# during the baseline period
use_yrs <- burn_in:baseline_yrs
temporal_model <- lm(temporal_sim$pcgr$Pop[use_yrs]~temporal_sim$temperature[use_yrs] + 
                       I(temporal_sim$temperature[use_yrs]^2))
# TO DO: replace this linear model with a Gompertz model



