function rho_m = rho_m2Layer(R_planet,dm,rho_av)

rho_w = 1000; % kg m-3

% Titan
LOCAL = 0;
if LOCAL ~= 0
	R_planet = 2575e3; % m
	dm = [50 200]*1e3
	rho_av = 2000; % kg m-3
end
rm = R_planet-dm; % m

rho_m = (R_planet.^3*rho_av-rho_w.*(R_planet.^3-rm.^3))./rm.^3;