function [P_bar,rho_rock] = d2P_2layer(R_planet_km,D_ocean_km,d_km,rho_av)
% Returns the pressure corresponding to the given depth in a hydrostatic ocean planet
% depth should be in m
% pressure is in Pa 
% from Turcotte and Schubert, page 86
rho_water = 1000; % kg/m3,  
R_planet_m = R_planet_km*1e3;
D_ocean_m = D_ocean_km*1e3;
d_m = d_km*1e3;

R_ocean_max_m = R_planet_m; % m
R_mantle_max_m = (R_planet_m-D_ocean_m); % m

if ~exist('rho_av')
	rho_av = 2000;
end
rho_rock = rho_m2Layer(R_planet_m,D_ocean_m,rho_av);
%rho_rock = 3800; % The Europa Model (Anderson et al, 1997) assuming no core
G = 6.67300e-11; % m3 kg-1 s-2

% P in the mantle
if find(d_m>D_ocean_m)
    dinds_mantle = find(d_m>=D_ocean_m);
    r_mantle_m = (R_planet_m-d_m(dinds_mantle));
    
    P_Pa(dinds_mantle) = 2/3*pi*G*(rho_rock^2*(R_mantle_max_m^2-r_mantle_m.^2)+...
                        rho_water^2*(R_ocean_max_m^2-R_mantle_max_m.^2)+...
                        2*rho_water*R_mantle_max_m^3*(rho_rock-rho_water)*(1/R_mantle_max_m-1/R_ocean_max_m));
end
% P in the ocean
if find(d_m<=D_ocean_m)
    test = find(d_m == D_ocean_m);
    if test
        d(test) = D_ocean_m - 1e-6;
    end
    dinds_ocean = find(d_m<=D_ocean_m);
    r_ocean_m = (R_planet_m-d_m(dinds_ocean));
   
    P_Pa(dinds_ocean) = 2/3*pi*G*(2*rho_water*R_mantle_max_m.^3.*(rho_rock-rho_water).*(1./r_ocean_m-1/R_ocean_max_m)...
                        +rho_water^2*(R_ocean_max_m^2-r_ocean_m.^2));
end

Pa2bar = 1e-5;
 P_bar = P_Pa*Pa2bar;


