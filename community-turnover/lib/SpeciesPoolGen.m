%
%
%	Help function to generate traits within the species pool
%	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu
%	Co-authors (and project leaders): J. Alexander - jake.alexander@unil.ch & L. Pellissier - loic.pellissier@usys.ethz.ch


function tr = SpeciesPoolGen(N,Tmin, Tmax, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax)

tr = zeros(N, 4) ;
%Competition intra
tr(:,1) = sort(rand(N,1)*(Tmax-Tmin) + Tmin);
tr(:,4) = rand(N,1)*(Cmax-Cmin) + Cmin;

%Growth rate
slope = (Gmax-Gmin)./(Tmax-Tmin);
tr(:,2) = tr(:,1).*slope+Gmin - slope*Tmin;

%Biomasse tolerance
slope = -(Lmax-Lmin)./(Tmax-Tmin);
tr(:,3) = tr(:,1).*slope+Lmax - slope*Tmin;

end


