#
#
# Help function to generate traits within the species pool
#	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu
#	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch

# Translated from Matlab to R by Peter Adler - peter.adler@usu.edu

SpeciesPoolGen <- function(N,Tmin, Tmax, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax){
  
  tr = matrix(0, N, 4) 
  
  # Competition intra
  tr[,1] <- sort(runif(N)*(Tmax-Tmin) + Tmin)
  tr[,4] <- runif(N)*(Cmax-Cmin) + Cmin
  
  # Growth rate
  slope <- (Gmax-Gmin)/(Tmax-Tmin)
  tr[,2] <- tr[,1]*slope+Gmin - slope*Tmin

  # Biomass tolerance
  slope <- -(Lmax-Lmin)/(Tmax-Tmin)
  tr[,3] <- tr[,1]*slope+Lmax - slope*Tmin

  return(tr)

}
  

