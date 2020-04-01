#
#
#	Help function to run community dynamics
#	Main script author: Loic Chalmandrier - lchalman@uwyo.edu, https://github.com/LoicChr
#	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch

# Translated from Matlab to R by Peter Adler - peter.adler@usu.edu

CommunityTempDis <- function(spxp_ini, tr, seed_rain,temp,dt){
#UNTITLED Summary of this function goes here
#   Detailed explanation goes here
  N <- dim(tr)[1]
  L_land <- NROW(spxp_ini)
  spxp_int <- seed_rain%*%spxp_ini
  Biom <- rowSums(spxp_int)
  
  growth_rate_temp <- (matrix(1,L_land,1)%*%t(tr[,2]))*(temp%*%matrix(1,1,N) - matrix(1,L_land,1)%*%t(tr[,1]))
  growth_rate_intra <- (matrix(1,L_land,1)%*%t(tr[,4]))*spxp_ini
  growth_rate_biom <- (matrix(1,L_land,1)%*%t(tr[,3]))*(Biom%*%matrix(1,1,N))
   
  spxp_final <- spxp_int + dt*spxp_int*(growth_rate_temp - growth_rate_intra - growth_rate_biom)
   
  spxp_final[spxp_final<10e-12] <- 0

  return(spxp_final)
  
}

generateTemps <- function(meanT,changeT,stationary=T){
  
  if(stationary==T){
    
    out <- meanT + c(rep(0,burnin_yrs+baseline_yrs),seq(changeT/warming_yrs,changeT,changeT/warming_yrs),
                        rep(changeT,final_yrs))
  }else{
    
    out <- meanT + seq(changeT/sim_yrs,changeT,changeT/sim_yrs)
    
  }
  
  return(out)
  
}
