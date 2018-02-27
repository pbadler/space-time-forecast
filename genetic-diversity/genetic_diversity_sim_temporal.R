
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
  N[iYr,] <- out$SB  # total seedbank
  Plants[iYr,] <- out$Plants  # germinated plants
  Fec[iYr,] <- out$SP # seed production
}   # next i


# Fit a "temporal" model describing abundance as a function of annual temperature
# during the baseline period. Assume Gompertz population growth.
N <-rowSums(N)
N_t0_baseline <- N[burn_in:(baseline_yrs-1)]
N_t1_baseline <- N[(burn_in+1):baseline_yrs]
r_baseline <- log(N_t1_baseline/N_t0_baseline)
temperature_t1_baseline <- temperature[(burn_in+1):baseline_yrs]

temporal_model <- lm(r_baseline ~ log(N_t0_baseline) + temperature_t1_baseline + 
                       I(temperature_t1_baseline^2))




