%
%
%	Example script to run the model of Alexander et al. 2018 'Lags in the response of mountain plant communities to climate change' - Global Change Biology
%	Main script author: Loïc Chalmandrier - lchalman@uwyo.edu
%


% Sourcing the help functions
addpath('./lib')

%generation of the Landscape aka a table with the following characteristics: X coordinate, Y coordinate, temperature
L_land = 500;
landscape = [1:L_land ; ones(1 ,L_land); linspace(0,20, L_land)]';

% Parameters to generate the species pool
N = 100;
Gmax = 0.5;
Gmin = 0.2;
Lmax = 1.5;
Lmin = 0.7;
Cmax = 0.2;
Cmin = 0.2; 
d = 10^-3;
 
% 1st step: initialisation of the metacommunity
   Tmax = 40000; dt = 0.025;
   spxp = ones(L_land, N).*(2/N);
   d = d./dt;

 % Generation of the species pool
 tr = SpeciesPoolGen(N, 4, 16.5, Gmin, Gmax, Lmin, Lmax, Cmin, Cmax);
 temp = landscape(:,3);
 
 % dispersal matrix
   seed_rain =  squareform(pdist(landscape(:,1:2)) < 2).*(d*dt*0.5);
   indx = logical( diag(ones(1,size(landscape,1)),0) );
   seed_rain(indx) = 1-d*dt;
   seed_rain(1,1) = 1-dt*d/2;
   seed_rain(size(landscape,1),size(landscape,1)) = 1-dt*d/2;
   
 %Run of the initialisation of community
 for t = 1:Tmax
     spxp = CommunityTempDis(spxp, tr, seed_rain, temp, dt);
 end 
 spxp_ini = spxp;

 %save the initial metacommunity
 mkdir('simulations/run_1')
 dlmwrite(strcat('simulations/run_1/spxp_ini.txt'),spxp_ini)
dlmwrite(strcat('simulations/run_1/landscape.txt'),landscape)
dlmwrite(strcat('simulations/run_1/tr.txt'),tr)

spxp_ini(spxp_ini < 10e-12) = 0;

spxp = spxp_ini;
 %From the initial situation, re-run community dynamic while temperature is increasing.
 temp = landscape(:,3);
 time_to_max = 20000;
 temp_increase = 3;
 for t = 1:time_to_max
     temp_new = temp_increase*t/time_to_max + temp;
     spxp = CommunityTempDis(spxp, tr, seed_rain, temp_new, dt);
 end 
 spxp_final = spxp;

% Save final meta-community
 dlmwrite(strcat('simulations/run_1/spxp_final.txt'),spxp_final)

