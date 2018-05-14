%
%
%	Help function to run community dynamics
%	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu
%	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch

function spxp_final = CommunityTempDis(spxp_ini, tr, seed_rain,temp,dt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    N = size(tr, 1);
  
   spxp_int = seed_rain*spxp_ini;
   Biom = sum(spxp_int,2);
   
   growth_rate_temp = (ones(500,1)*tr(:,2)').*(temp*ones(1,N) -ones(500,1)*tr(:,1)');
   growth_rate_intra = (ones(500,1)*tr(:,4)').*spxp_ini;
   growth_rate_biom = (ones(500,1)*tr(:,3)').*(Biom*ones(1,N));
   
   spxp_final = spxp_int + dt.*spxp_int.*(growth_rate_temp - growth_rate_intra - growth_rate_biom);
   
   spxp_final(spxp_final<10e-12) = 0;

end

