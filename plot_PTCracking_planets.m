 function Planets = plot_PTCracking_planets
clear
global R_planet_m d_ocean_m d_ocean_vector_m Tdot_str...
  n Tprime To IS_EARTH IS_MARS USE_P_EFFECTIVE rho_av rho_core_max... 
     PLOT_ZT_1MM PLOT_ZT_10MM PLOT_HT PLOT_H2 PLOT_GRL PLOT_AB...
	t_style d_Earth P_PREM P_Earth T_Earth POLY_ORDER lambda PLOT_ALL_PT...
    PLOT_ERRZ_VS_R INCLUDE_LARGER_ICY_WORLDS PLOT_SCHEMATIC

PROCESS_NEW = 1;


INCLUDE_MAGENTAS = 1;
INCLUDE_LARGER_ICY_WORLDS = 1;
if INCLUDE_LARGER_ICY_WORLDS
    larger_icy_planets_plot_color = 'y';
else
    larger_icy_planets_plot_color = 'n';
end

% brittle ductile criteria taken from Fournier 1996.  These would be higher
% for basalt according to the author (though no specific values are supplied), so the values assumed are conservative
P_bd_MPa = 120;
T_bd_oC = 325; 

USE_P_EFFECTIVE = 0;
PLOT_ALL_PT = 1;
PLOT_ZT_1MM = 1; % 
PLOT_ZT_10MM = 0;
PLOT_ERRZ_VS_R = 0;
PLOT_HT = 1;
PLOT_H2 = 1; % if calculating H2 vs time, calculate Hserp and H2 generation and output tex table for AB paper
             % if not plotting H2, calculate cracking for current-era heating rates and plotting P vs T, output tex table for GRL paper

 PLOT_SCHEMATIC = 0;
             
if ~PLOT_H2
	t_solar_system = 0;
else
%     t_solar_system = logspace(0,9.648,200);
	t_solar_system=  [0 10.^(1:8) 10.^(8.5:0.5:9) (1.1:0.1:4.5)*1e9];
	%t_solar_system=  [0 10.^(0.1:0.01:9) (4.5)*1e9];
end

n = 1;

%     lamda = 2.0; % W m^-1 k^-1 Lowell and Dubose, 2005
%       lamda = 3.0; % W m^-1 k^-1 Hofmeister 1999

% for lambda = 3 % W / m / K
lmatrix = 3;
for hi = 1:length(lmatrix) %5.2 % W / m / K
    lambda = lmatrix(hi);
Tdot_str = '1oCyr';
%Tdot_str = '1oCMyr'; % yields a factor of 1/2 of 1oCyr
%Tdot_str = '1oCGyr'; % yields a factor of 1/4 of 1oCyr

outer_planets_plot_color = [100 0 255]/255; % RGB purple
To_ammonia = -20; % 10% ammonia, allowing for modest heat flow

if PROCESS_NEW 
deMartinCrackingDepthVGrainSize;

%H = H_U + H_Th + H_K = sum( delta_C * dE/dC )
%going back in time, concentration increases as C(t) = Co exp(lambda*t)
%where lambda_U = 1.551e-10 per yr, lambda_Th = 0.495e-10 per yr, and
%lambda_K = 5.543e-10 per yr
% % 
%============
t_planet = 'Earth';
IS_EARTH = 1;
R_planet_m = 6371e3; A_planet = 4*pi*R_planet_m^2;
d_ocean_vector_m = 4e3;
To = 0;
rho_av = 5520;
f_Si_Earth = 1;
M_Earth = 5.9742e24;% kg
% H = 1e13; %5e-3*A_planet; % Lorenz 2002, page 179, last paragraph; source not cited
H_Earth = get_pastRadiogenicHeat_UThK(M_Earth,f_Si_Earth,t_solar_system);
Earth = find_PT(t_planet,'g',R_planet_m,M_Earth,rho_av,d_ocean_vector_m,H_Earth,To,t_solar_system)
Planets(n) = Earth;
IS_EARTH = 0;
% % 
n = n+1;
%============
t_planet = 'Mars';
IS_MARS = 1;
R_planet_m = 3397e3; A_planet = 4*pi*R_planet_m^2;
d_ocean_vector_m = [0]*1e3;
To = 0;
rho_av = 3930;
M_Mars = 6.42e23; % kg
f_Si_Mars = 1;
t_Mars = t_solar_system;
H_Mars = get_pastRadiogenicHeat_UThK(M_Mars,f_Si_Mars,t_Mars);
Mars =find_PT(t_planet, 'r',R_planet_m,M_Mars,rho_av,d_ocean_vector_m,H_Mars,To,t_solar_system)
Mars.t_recharge_Mya = 3.5e3;
Planets(n) = Mars;
IS_MARS = 0;

% Science 11 April 2003:
% Vol. 300. no. 5617, pp. 299 - 303
% DOI: 10.1126/science.1079645
% 	
% Prev | Table of Contents | Next
% Research Articles
% 
% Fluid Core Size of Mars from Detection of the Solar Tide
% 
% C. F. Yoder,* A. S. Konopliv,* D. N. Yuan, E. M. Standish, W. M. Folkner
% 
% The solar tidal deformation of Mars, measured by its k2 potential Love number, has been obtained from an analysis of Mars Global Surveyor radio tracking. 
% The observed k2 of 0.153 ± 0.017 is large enough to rule out a solid iron core and so indicates that at least the outer part of the core is liquid. 
% The inferred core radius is between 1520 and 1840 kilometers and is independent of many interior properties, although partial melt of the mantle is one factor that could reduce core size. 
% Ice-cap mass changes can be deduced from the seasonal variations in air pressure and the odd gravity harmonic J3, given knowledge of cap mass distribution with latitude. The south cap seasonal mass change is about 30 to 40% larger than that of the north cap.

n = n+1;
%==========
% info from Wikipedia
t_planet = 'Ceres';
R_planet_m = 475e3; % pm 3e3 m 
%  d_ocean_vector_m = [52]*1e3; % m
d_ocean_vector_m = [100]*1e3; % m
To = 273.15-273.15;
rho_av = 2206;
M_Ceres = 9.46e20; %  ± 0.04  kg
% f_Si_Ceres = 0.94;% to 1; Schubert et al 1986 for Europa
f_Si_Ceres = get_SiMassFraction(rho_av,1000,2900); % 
t_Ceres = t_solar_system;
H_Ceres = get_pastRadiogenicHeat_UThK(M_Ceres,f_Si_Ceres,t_Ceres);
Ceres = find_PT(t_planet,'b',R_planet_m,M_Ceres,rho_av,d_ocean_vector_m,H_Ceres,To,t_solar_system)
Planets(n) = Ceres;
% 
% n = n+1;
%  d_ocean_vector_m = [120]*1e3; % m
% t_planet = 'Ceres 120';
% Ceres = find_PT(t_planet,'b',R_planet_m,M_Ceres,rho_av,d_ocean_vector_m,H_Ceres,To,t_solar_system)
% Planets(n) = Ceres;

 n = n+1;
%==========
R_planet_m = 1565e3; % m
% d_ocean_vector_m = [80 170]*1e3; % m
d_ocean_vector_m = [100]*1e3; % m
To = 0;
rho_av = 2970;
M_Europa = 4870e19; % kg
f_Si_Europa = 0.94;% to 1; Schubert et al 1986
t_Europa = t_solar_system;
t_planet = 'Europa';
H_Europa = get_pastRadiogenicHeat_UThK(M_Europa,f_Si_Europa,t_Europa);
Europa = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa,To,t_solar_system)
Europa.t_recharge_Mya = 2e3;
Planets(n) = Europa;

% n= n+1;
% t_planet = 'Europa 170';
% d_ocean_vector_m = [170]*1e3; % m
% H_Europa = get_pastRadiogenicHeat_UThK(M_Europa,f_Si_Europa,t_Europa);
% Europa = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa,To,t_solar_system)
% Planets(n) = Europa;
% t_planet = 'Europa (radiogenic)';
% H_Europa = get_pastRadiogenicHeat_UThK(M_Europa,f_Si_Europa,t_Europa);
% Europa = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa,To,t_solar_system)
% Planets(n) = Europa;
% 
% n = n+1; 
%  t_planet = 'Europa (intermediate_0)';
% %t_planet = '';
% H_tidal  = 0.35e12; % from Vance et al. 2007
% H_Europa0 = H_Europa + H_tidal;
% Europa0 = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa0,To,t_solar_system)
% Planets(n) = Europa0;
% 
% 
% n = n+1; 
% % t_planet = 'Europa (intermediate_1)';
% t_planet = '';
% H_tidal  = 0.7e12; % from Vance et al. 2007
% H_Europa10 = H_Europa + H_tidal;
% Europa10 = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa10,To,t_solar_system)
% Planets(n) = Europa10;
% 
% n = n+1; 
% % t_planet = 'Europa (intermediate_2)';
% t_planet = '';
% H_tidal  = 1.4e12; % from Vance et al. 2007
% H_Europa15 = H_Europa + H_tidal;
% Europa15 = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa15,To,t_solar_system)
% Planets(n) = Europa15;
% 
% n = n+1; 
% t_planet = 'Europa (tidal)';
% H_tidal  = 2.8e12; % from Vance et al. 2007
% H_Europa2 = H_Europa + H_tidal;
% Europa2 = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa2,To,t_solar_system)
% Europa2.t_recharge_Mya = 3e3;
% Planets(n) = Europa2;

% n = n+1; 
% t_planet = 'Europa (tidal_2)';
% H_tidal  = 3.5e12; % Lowell and Dubose 2005
% H_Europa2 = H_Europa + H_tidal;
% Europa3 = find_PT(t_planet,'b',R_planet_m,M_Europa,rho_av,d_ocean_vector_m,H_Europa2,To,t_solar_system)
% Europa3.t_recharge_Mya = 3e3;
% Planets(n) = Europa3;

%hydrothermal fluid flux from Lowell and Dubose (2005) is estimated at 10^7
%kg/s. Average black smoker H2 flux is 7 mmol/kg (Holland et al. 2002)
% 7e4*pi*1e7
%      2.199114857512855e+12 mol/yr globally
%                         =    4.4150e+18 molecules/cm2/yr

if INCLUDE_LARGER_ICY_WORLDS
n = n+1;
%===========
t_planet = 'Callisto';
R_planet_m = 2403e3;
d_ocean_vector_m = 450*1e3;
To = 0;
rho_av = 1850;
M_Callisto =  1.35e23; % kg
t_Callisto = t_solar_system;
f_Si_Callisto = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Callisto = get_pastRadiogenicHeat_UThK(M_Callisto,f_Si_Callisto,t_Callisto);
Callisto = find_PT(t_planet, larger_icy_planets_plot_color,R_planet_m,M_Callisto,rho_av,d_ocean_vector_m,H_Callisto,To,t_solar_system)
Planets(n) = Callisto;
% figure(370);plot_tz(t_solar_system/1e6,Callisto.z_cracking_1mm_m/1e3,t_planet,outer_planets_plot_color);
% figure(371);plot_tz(t_solar_system/1e6,Callisto.z_cracking_10mm_m/1e3,t_planet,outer_planets_plot_color);
% 

n = n+1;
%===========
t_planet = 'Ganymede';
R_planet_m = 2631.2e3;
d_ocean_vector_m = 900*1e3;
To = 50;
rho_av = 1942;
M_Ganymede =  1.35e23; % kg
t_Ganymede = t_solar_system;
f_Si_Ganymede = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Ganymede = get_pastRadiogenicHeat_UThK(M_Ganymede,f_Si_Ganymede,t_Ganymede);
Ganymede = find_PT(t_planet, larger_icy_planets_plot_color,R_planet_m,M_Ganymede,rho_av,d_ocean_vector_m,H_Ganymede,To,t_solar_system)
Planets(n) = Ganymede;


n = n+1;
%===========
t_planet = 'Kepler 78b';
R_planet_m = 1*6371e3; %1.5 times the radius of Earth, referring to Zeng et al. 2014
d_ocean_vector_m = 800*1e3;
To = 0;
rho_av = 1942;
M_Kepler78b =  1*M_Earth; % kg
t_Kepler78b = t_solar_system;
f_Si_Kepler78b = get_SiMassFraction(rho_av,1400,5000); 
H_Kepler78b = get_pastRadiogenicHeat_UThK(M_Kepler78b,f_Si_Kepler78b,t_Kepler78b);
Kepler78b = find_PT(t_planet, larger_icy_planets_plot_color,R_planet_m,M_Kepler78b,rho_av,d_ocean_vector_m,H_Kepler78b,To,t_solar_system)
Planets(n) = Kepler78b;
end


n = n+1;
%==========
t_planet = 'Enceladus';
R_planet_m = 252e3; 
d_ocean_vector_m = 80*1e3;
To = To_ammonia;
% densities from Porco et al 2006
rho_av = 1608;
f_Si_Enceladus = 0.6; % rho_Si = 2500-3500 gives f_Si = 0.5617=0.6600
M_Enceladus = 7.30e19;
t_Enceladus = t_solar_system;
H_Enceladus = get_pastRadiogenicHeat_UThK(M_Enceladus,f_Si_Enceladus,t_Enceladus);
Enceladus = find_PT(t_planet, 'c',R_planet_m,M_Enceladus,rho_av,d_ocean_vector_m,H_Enceladus,To,t_solar_system)
Enceladus.z_cracking_1mm_m = ones(1,length(t_solar_system))*(R_planet_m-d_ocean_vector_m);
Enceladus.z_cracking_10mm_m = ones(1,length(t_solar_system))*(R_planet_m-d_ocean_vector_m);
Enceladus.P_z_1mm_MPa = Enceladus.P_MPa(end);Enceladus.P_z_10mm_MPa = Enceladus.P_MPa(end);
Enceladus.d_ocean_km = d_ocean_vector_m/1e3;
Planets(n) = Enceladus;

%  n = n+1;
% % %==========
% t_planet = 'Dione';
% R_planet_m = 562.5e3; 
% d_ocean_vector_m = Inf;
% To = To_ammonia;%   % densities from Porco et al 2006
% rho_av = 1470;
% f_Si_Dione = 0.6; % rho_Si = 2500-3500 gives f_Si = 0.5617=0.6600
% M_Dione = 10.96e20;
% t_Dione = t_solar_system;
% H_Dione = get_pastRadiogenicHeat_UThK(M_Dione,f_Si_Dione,t_Dione);
% Dione = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Dione,rho_av,d_ocean_vector_m,H_Dione,To,t_solar_system)
% Planets(n) = Dione;

 n = n+1;
%==========
t_planet = 'Rhea';
R_planet_m = 764.5e3; 
d_ocean_vector_m = 417.3e3;
To = To_ammonia;%176-273.15;
% densities from Porco et al 2006
rho_av = 1234;
f_Si_Rhea = 0.6; % rho_Si = 2500-3500 gives f_Si = 0.5617=0.6600
M_Rhea = 23.10e20;
t_Rhea = t_solar_system;
H_Rhea = get_pastRadiogenicHeat_UThK(M_Rhea,f_Si_Rhea,t_Rhea);
Rhea = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Rhea,rho_av,d_ocean_vector_m,H_Rhea,To,t_solar_system)
Planets(n) = Rhea;

n = n+1;
%==========
t_planet = 'Iapetus';
R_planet_m = 736e3; % mean from wiki overall shape section
d_ocean_vector_m = 100e3;
To = To_ammonia;%176-273.15;
% densities from Porco et al 2006
rho_av = 1083;
f_Si_Iapetus = 0.6; % rho_Si = 2500-3500 gives f_Si = 0.5617=0.6600
M_Iapetus = 1.8e21;
t_Iapetus = t_solar_system;
H_Iapetus = get_pastRadiogenicHeat_UThK(M_Iapetus,f_Si_Iapetus,t_Iapetus);
Iapetus = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Iapetus,rho_av,d_ocean_vector_m,H_Iapetus,To,t_solar_system)
Iapetus.z_cracking_1mm_m = ones(1,length(t_solar_system))*(R_planet_m-d_ocean_vector_m);
Planets(n) = Iapetus;

if INCLUDE_LARGER_ICY_WORLDS
n = n+1;
%===========
t_planet = 'Titan';
R_planet_m = 2575e3;
d_ocean_vector_m = [600]*1e3;
To = To_ammonia;
rho_av = 1880;
M_Titan =  1.3455e23; % kg
t_Titan = t_solar_system;
f_Si_Titan = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Titan = get_pastRadiogenicHeat_UThK(M_Titan,f_Si_Titan,t_Titan);
Titan = find_PT(t_planet, larger_icy_planets_plot_color,R_planet_m,M_Titan,rho_av,d_ocean_vector_m,H_Titan,To,t_solar_system)
Planets(n) = Titan;
end

if INCLUDE_MAGENTAS
n = n+1;
%===========
t_planet = 'Miranda';
R_planet_m = 235.8e3; %±.7km Thomas, P. C. (1988). "Radii, shapes, and topography of the satellites of Uranus from limb coordinates". Icarus 73 (3): 427?441. Bibcode:1988Icar...73..427T. doi:10.1016/0019-1035(88)90054-1. edit 
d_ocean_vector_m = Inf;
To = To_ammonia;
rho_av = 1200;
M_Miranda =  6.59e19; % kg Jacobson, R. A.; Campbell, J. K.; Taylor, A. H.; Synnott, S. P. (June 1992). "The masses of Uranus and its major satellites from Voyager tracking data and earth-based Uranian satellite data". The Astronomical Journal 103 (6): 2068?2078. Bibcode:1992AJ....103.2068J. doi:10.1086/116211
f_Si_Miranda = get_SiMassFraction(rho_av,1000,2500); % 
H_Miranda = get_pastRadiogenicHeat_UThK(M_Miranda,f_Si_Miranda,t_solar_system);
Miranda = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Miranda,rho_av,d_ocean_vector_m,H_Miranda,To,t_solar_system)
Planets(n) = Miranda;    
    
n = n+1;
%===========
t_planet = 'Ariel';
R_planet_m = 578.9e3; 
d_ocean_vector_m = Inf;
To = To_ammonia;
rho_av = 1665;
M_Ariel =  13.53e20; % kg 
f_Si_Ariel = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice --> this is probably a copy and paste error from the Titan or other case
H_Ariel = get_pastRadiogenicHeat_UThK(M_Ariel,f_Si_Ariel,t_solar_system);
Ariel = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Ariel,rho_av,d_ocean_vector_m,H_Ariel,To,t_solar_system)
Planets(n) = Ariel;

n = n+1;
%===========
t_planet = 'Umbriel';
R_planet_m = 584.7e3; 
d_ocean_vector_m = Inf;
To = To_ammonia;
rho_av = 1400;
M_Umbriel =  11.72e20; % kg
f_Si_Umbriel = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Umbriel = get_pastRadiogenicHeat_UThK(M_Umbriel,f_Si_Umbriel,t_solar_system);
Umbriel = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Umbriel,rho_av,d_ocean_vector_m,H_Umbriel,To,t_solar_system)
Planets(n) = Umbriel;

n = n+1;
%===========
t_planet = 'Titania';
R_planet_m = 788.9e3; 
d_ocean_vector_m = 269.2e3;
To = To_ammonia;%204-273.15;
rho_av = 1715;
M_Titania =  35.27e20; % kg
f_Si_Titania = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Titania = get_pastRadiogenicHeat_UThK(M_Titania,f_Si_Titania,t_solar_system);
Titania = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Titania,rho_av,d_ocean_vector_m,H_Titania,To,t_solar_system)
Planets(n) = Titania;

n = n+1;
%===========
t_planet = 'Oberon';
R_planet_m = 761.4e3; 
d_ocean_vector_m = 280.3e3;
To = To_ammonia;%194-273.15;
rho_av = 1630;
M_Oberon =  30.14e20; % kg
f_Si_Oberon = get_SiMassFraction(rho_av,1200,2900); % use high pressure density of ice
H_Oberon = get_pastRadiogenicHeat_UThK(M_Oberon,f_Si_Oberon,t_solar_system);
Oberon = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Oberon,rho_av,d_ocean_vector_m,H_Oberon,To,t_solar_system)
Planets(n) = Oberon;
end

if INCLUDE_LARGER_ICY_WORLDS
n = n+1;
%===========
t_planet = 'Triton';
R_planet_m = 1353.4e3; 
d_ocean_vector_m = 336.4e3;
To = To_ammonia; %255-273.15;
rho_av = 2061;
M_Triton =  214.0e20; % kg
f_Si_Triton = 0.65; % use high pressure density of ice
H_Triton = get_pastRadiogenicHeat_UThK(M_Triton,f_Si_Triton,t_solar_system);
Triton = find_PT(t_planet, larger_icy_planets_plot_color,R_planet_m,M_Triton,rho_av,d_ocean_vector_m,H_Triton,To,t_solar_system)
Planets(n) = Triton;
end

if INCLUDE_MAGENTAS
% 
n = n+1;
%===========
t_planet = 'Pluto';
R_planet_m = 1195.0e3;% 1189
d_ocean_vector_m = 364.8e3;
To = To_ammonia; %252-273.15;
rho_av = 1838; % 1863
M_Pluto =  131.4e20; % kg
f_Si_Pluto = 0.65; % use high pressure density of ice
H_Pluto = get_pastRadiogenicHeat_UThK(M_Pluto,f_Si_Pluto,t_solar_system);
Pluto = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Pluto,rho_av,d_ocean_vector_m,H_Pluto,To,t_solar_system)
Planets(n) = Pluto;

n = n+1;
%===========
t_planet = 'Charon';
R_planet_m = 603.6e3; 
d_ocean_vector_m = Inf;
To = To_ammonia;
rho_av = 1757; % 1707
M_Charon =  16.2e20; % kg
f_Si_Charon = 0.65; % use high pressure density of ice
H_Charon = get_pastRadiogenicHeat_UThK(M_Charon,f_Si_Charon,t_solar_system);
Charon = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Charon,rho_av,d_ocean_vector_m,H_Charon,To,t_solar_system)
Planets(n) = Charon;

n = n+1;
%===========
t_planet = 'Eris'; % Formerly '2003 UB_{313}'
R_planet_m = 1300e3; 
d_ocean_vector_m = 386e3;
To = To_ammonia;
rho_av = 1757;
M_Eris =  172e20; % kg
f_Si_Eris = 0.65; % use high pressure density of ice
H_Eris = get_pastRadiogenicHeat_UThK(M_Eris,f_Si_Eris,t_solar_system);
Eris = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Eris,rho_av,d_ocean_vector_m,H_Eris,To,t_solar_system)
Planets(n) = Eris;

n = n+1;
%===========
t_planet = 'Sedna / Orcus'; % Orcus was formerly 2004 DW
R_planet_m = 800e3; 
d_ocean_vector_m = 238e3;
To = To_ammonia;
rho_av = 1757;
M_Sedna_Orcus =  40e20; % kg
f_Si_Sedna_Orcus = 0.65; % Hussmann2006 page 12, same as for Triton, Pluto and Charon
H_Sedna_Orcus = get_pastRadiogenicHeat_UThK(M_Sedna_Orcus,f_Si_Sedna_Orcus,t_solar_system);
Sedna_Orcus = find_PT(t_planet, outer_planets_plot_color,R_planet_m,M_Sedna_Orcus,rho_av,d_ocean_vector_m,H_Sedna_Orcus,To,t_solar_system)
Planets(n) = Sedna_Orcus;
end

save('Planets.mat','Planets')
else
    Planets = load('Planets.mat');
    Planets = Planets.Planets;
end
Planets = find_bd(Planets);

if PLOT_H2
    Planets = SerpHeatOceanPlanets(Planets);
%     for ij = 1:length(Planets)
%         name{hi,ij} = Planets(ij).name;
%         if ~isempty(Planets(ij).F_serp_W_m2)
%             Fserp(hi,ij) = Planets(ij).F_serp_W_m2(1);
%             FH2(hi,ij)  = Planets(ij).FH2_molecules_cm2_s(1);
%         end            
%             meanFserp(hi,ij) = Planets(ij).mean_F_serp_W_m2;
%             meanFH2(hi,ij) = Planets(ij).mean_FH2_molecules_cm2_s;
%     end
% 
% disp([Fserp' meanFserp']);

if PLOT_ZT_1MM
figure(370);clf;set(gcf,'Name','z vs t, 1mm');
xlabel('Time Before Present (Ga)','FontSize',36,'FontWeight','bold');
ylabel('Cracking Depth (km)','FontSize',36,'FontWeight','bold');
set(gca,'XDir','reverse','YAxisLocation','right','YScale','log','YLim',[0.6 250],'YGrid','on','FontSize',20);
box on;hold on;

figure(372);set(gcf,'Name','dzdt vs t, 1mm');
xlabel('Time Before Present (Mya)','FontSize',36,'FontWeight','bold');
ylabel('dzdt (mm yr^{-1})','FontSize',36,'FontWeight','bold');
set(gca,'XDir','reverse','YAxisLocation','right','YScale','log','YGrid','on','FontSize',20);
box on;hold on;

figure(373);clf;set(gcf,'Name','z vs R, 1mm');
xlabel('R_{planet} (km)','FontSize',36,'FontWeight','bold');
ylabel('Cracking Depth (km)','FontSize',36,'FontWeight','bold');
set(gca,'YAxisLocation','right','YScale','log','YGrid','on','FontSize',20);
box on; hold on;
end

if PLOT_ZT_10MM
figure(371);clf;set(gcf,'Name','z vs t, 10mm');
xlabel('Time Before Present (Ga)','FontSize',36,'FontWeight','bold');
ylabel('Cracking Depth (km)','FontSize',36,'FontWeight','bold');
set(gca,'XDir','reverse','YAxisLocation','right','YScale','log','YGrid','on','FontSize',20);
box on; hold on;

figure(374);clf;set(gcf,'Name','z vs R, 10mm');
xlabel('R_{planet} (km)','FontSize',36,'FontWeight','bold');
ylabel('Cracking Depth (km)','FontSize',36,'FontWeight','bold');
set(gca,'YAxisLocation','right','YScale','log','YGrid','on','FontSize',20);
box on; hold on;
end

if PLOT_HT
figure(372);set(gcf,'Name','H vs t');clf;
xlabel('Time Before Present (Ga)','Fontsize',36,'FontWeight','bold');
% ylabel('Internal Radiogenic Heating (W)','FontSize',36,'FontWeight','bold');
set(gca,'XDir','reverse','YScale','log','YAxisLocation','right','YGrid','on','FontSize',20);
box on; hold on;
Planets = plot_tz(Planets);
end

end
figure(369);set(gcf,'Name','P vs T, 1mm'); 
set(gca,'YDir','reverse','YLim',[0 Inf],'XLim',[-90 550],'Box','on','FontSize',20);
title('Thermal Cracking','FontSize',36,'FontWeight','bold')
xlabel('T (^oC)','FontSize',36,'FontWeight','bold');
ylabel('P_c (MPa)','FontSize',36,'FontWeight','bold');
plot_PT(Planets);
 save_to_table(Planets);
    if PLOT_SCHEMATIC
        plot_schematics(Planets)
    end


end
 end
%==========================================================================
function Planet_struct = find_PT(t_name,t_color,R_planet_m,M_planet,rho_av,d_ocean_vector_m,H,To,time)
global CrackingFront IS_EARTH P_PREM P_Earth
clear rho_mc;
rho_mc = ones(2,length(d_ocean_vector_m));

if d_ocean_vector_m == Inf; % if C/MR2 data are not available to constrain the "ocean" depth, make a conservative estimate based on a "mantle" density of 3000 kg m-3
	rho_mc = 2900;
	rho_w = 1000;
	d_ocean_vector_m = R_planet_m*(1-((rho_av-rho_w)/(rho_mc-rho_w))^(1/3));
end
	for ij = 1:length(d_ocean_vector_m)
		for jk = 1:length(H) % for every H (i.e., for heat H determined at time "time")
			d_ocean_m = d_ocean_vector_m(ij);
			[P(ij,:),T(ij,:,jk),rho_mc(ij),z_m] = get_PT(d_ocean_m,R_planet_m,rho_av,H(jk),To);
        end
%         for jk = 1:length(squeeze(T(ij,:,1))) % for every depth below the seafloor
%             pp = spline(time,squeeze(T(ij,jk,:)));
%             Tdot(ij,jk,:) = pp1derv(pp,time); % calculate the cooling rate in oC per year
%         end
	end

d_inds = find(rho_mc<6000 & rho_mc>2700);

%Density of Olivine is (av) 3250, Density of Serpentine is 
if ~isempty(d_inds)
	d_ocean_m = d_ocean_vector_m(d_inds); % m
	rho_mc = rho_mc(d_inds); P = P(d_inds,:); T = T(d_inds,:,:); 

    for ij = 1:length(d_ocean_m)
        max_z_m = R_planet_m-d_ocean_m(ij);
        P_planet_MPa = P(ij,:);
        for jk = 1:length(H) 
            T_planet_oC = squeeze(T(ij,:,jk));
            %This is the function:  
                %[P_cracking_MPa,T_cracking_oC,z_cracking_m] =
                %get_Pz_cracking(P_front_MPa,T_front,P_planet_MPa,T_planet,max_P_MPa,R_planet_m,d_ocean_m,max_z,rho_av)
			[P_cracking_p1mm(ij,jk),T_cracking_p1mm(ij,jk),z_cracking_p1mm(ij,jk)]=...
                get_Pz_cracking(CrackingFront.Pc_p1mm_MPa,...
                                CrackingFront.T_p1mm_oC,...
                                P_planet_MPa,T_planet_oC,...
                                max(CrackingFront.Pc_p1mm_MPa),...
                                R_planet_m,d_ocean_m(ij),max_z_m,rho_av);
            [P_cracking_1mm(ij,jk),T_cracking_1mm(ij,jk),z_cracking_1mm(ij,jk)]=...
                get_Pz_cracking(CrackingFront.Pc_1mm_MPa,...
                                CrackingFront.T_1mm_oC,...
                                P_planet_MPa,T_planet_oC,...
                                max(CrackingFront.Pc_1mm_MPa),...
                                R_planet_m,d_ocean_m(ij),max_z_m,rho_av);
			[P_cracking_10mm(ij,jk),T_cracking_10mm(ij,jk),z_cracking_10mm(ij,jk)]=...
                get_Pz_cracking(CrackingFront.Pc_10mm_MPa,...
                                CrackingFront.T_10mm_oC,...
                                P_planet_MPa,T_planet_oC,...e
                                max(CrackingFront.Pc_10mm_MPa),...
                                R_planet_m,d_ocean_m(ij),max_z_m,rho_av);
        end
    end
  
else
	z_cracking_p1mm = 0; T_cracking_p1mm = 0; P_cracking_p1mm = 0;
    z_cracking_1mm = 0; T_cracking_1mm = 0;  P_cracking_1mm = 0;
    z_cracking_10mm = 0; T_cracking_10mm = 0;  P_cracking_10mm = 0;
end

Planet_struct = struct('name',t_name,'plot_color',t_color,...
                        'z_cracking_p1mm_m',z_cracking_p1mm,'P_z_p1mm_MPa',P_cracking_p1mm,'T_z_p1mm_oC',T_cracking_p1mm,...
						'z_cracking_1mm_m',z_cracking_1mm,'P_z_1mm_MPa',P_cracking_1mm,'T_z_1mm_oC',T_cracking_1mm,...
						'z_cracking_10mm_m',z_cracking_10mm,'P_z_10mm_MPa',P_cracking_10mm,'T_z_10mm_oC',T_cracking_10mm,...
						'R_m',R_planet_m,'M_kg',M_planet,'d_ocean_km',d_ocean_m/1e3,'rho_av',rho_av,'rho_mc',rho_mc,...
                        'depth_km',z_m*1e-3, 'P_MPa',P,'T_oC',T,'H_W',H,'t_yr',time,'t_recharge_Mya',[]);% 'Tdot_oC_yr',Tdot,
end
%===============================================================================        
function Planets = find_bd(Planets)
    
% brittle ductile criteria taken from Fournier 1999.  These would be higher
% for basalt according to the author (though no specific values are supplied), so the values assumed are conservative
P_bd_MPa = 120; % sigma1-sigma3 equals this value at the brittle-ductile transition
T_bd_oC = 350; % the geothermal temperature equals or exceeds this value at the brittle-ductile transition
for iPlanet = 1:length(Planets)
        Planet = Planets(iPlanet);
        if P_bd_MPa>max(Planet.P_MPa)
            Planets(iPlanet).z_bd_m = [];
        else
       z_bd_m = spline(Planet.P_MPa,Planet.depth_km*1e3,P_bd_MPa);
       Planets(iPlanet).z_bd_m = ones(size(Planet.t_yr));
       Planets(iPlanet).Tz_bd_oC = ones(size(Planet.t_yr));
       for id = 1:length(Planet.T_oC(:,1,1)) % for each putative ocean depth
           for iT = 1:length(Planet.T_oC(1,1,:)) % for each time step
            T_oC = squeeze(Planet.T_oC(id,:,iT));
                Tz_bd_oC = spline(Planet.depth_km,T_oC,z_bd_m*1e-3);   
                    if Tz_bd_oC<T_bd_oC
                        if T_bd_oC > max(T_oC)
                            Planets(iPlanet).z_bd_m(id,iT) = NaN;
                        else
                           Tinds = T_oC>Tz_bd_oC+1;
                            Planets(iPlanet).z_bd_m(id,iT) = spline(T_oC(Tinds),Planet.depth_km(Tinds)*1e3,T_bd_oC);
                        end
                    else
                        Planets(iPlanet).z_bd_m(id,iT)  = z_bd_m;
                        Planets(iPlanet).Tz_bd_oC(id,iT)  = z_bd_m;
                    end
            end
       end
    end
end
end
%===============================================================================        
function plot_PT(Planets)
global PLOT_H2 PLOT_ALL_PT PLOT_ZT_1MM PLOT_ZT_10MM INCLUDE_LARGER_ICY_WORLDS
    figure(369);hold on;
    lP = length(Planets);
    n = 0;
    
    
    
    if INCLUDE_LARGER_ICY_WORLDS
        axis([-100 700 0 10000]);
    end
    
    for ij = 1:lP
        if PLOT_ALL_PT% & Planets(ij).z_cracking_1mm_m>0
            for jk = 1:length(Planets(ij).d_ocean_km)
                mantle_inds = find(Planets(ij).depth_km>=Planets(ij).d_ocean_km(jk));
%                 Planet_inds(n) = ij;
                T = Planets(ij).T_oC(jk,mantle_inds);
                P = Planets(ij).P_MPa(jk,mantle_inds);
                plot(T,P,'Color',Planets(ij).plot_color);%,[LineStyle{n}],'LineWidth',LineWidth)
        %    H = Planets(ij).H;
        % 	lH = log10(H);powerH = floor(lH);magH = 10^(lH-powerH); % format H to display as magH x 10^powerH
        % 	if time == 0;
        % 		t_heat = ['; H = ' num2str(magH,'%0.2f') ' \times 10^{' num2str(powerH,'%0.0f') '} W'];
        % 	else
        % 		t_heat = ['; t = ',num2str(time/1e6),' Mya; H = ' num2str(magH,'%0.2f') ' \times 10^{' num2str(powerH,'%0.0f') '} W'];
        % 	end
           d_text{1} = Planets(ij).name;
            if PLOT_ZT_1MM
                d_text = {d_text{:} ['     z_{1 mm} = '...
                            num2str(Planets(ij).z_cracking_1mm_m(jk)/1e3,'%0.0f') ' km']};
            end
            if PLOT_ZT_10MM %&& Planets(ij).z_cracking_10mm_m>0
%                 if Planets(ij).z_cracking_10mm_m<0 % set to zero if
%                 z_cracking is nonsense
%                     Planets(ij).z_cracking_10mm_m = Planets(ij).R_m-Planets(ij).d_ocean_km*1e3;
%                 end
                d_text = {d_text{:} ['     z_{10 mm} = ' ...
                    num2str(Planets(ij).z_cracking_10mm_m(jk)/1e3,'%0.0f') ' km']};
            end
            t= text(T(1),P(1),d_text,'FontSize',16);
            clear d_text
%             text_rot(t,T,P);
        end
    end
    end
end
 %=========================================================================
function plot_schematics(Planets)
    figure(368);set(gcf,'Name','Terrestrial Planets Schematic'); clf;
%     set(gca,'XLim',[0 20],'Box','on','FontSize',20);
    set(gca,'XLim',[0 252],'Box','on','FontSize',20);
 
    XLim = [-100 300];
    Pl_inds = [0 0 0 0];
    for iP = 1:length(Planets)
        if strcmp(Planets(iP).name,'Earth') 
            Pl_inds(1) = iP;
        elseif strcmp(Planets(iP).name,'Mars') 
            Pl_inds(2) = iP;
        elseif strcmp(Planets(iP).name,'Europa') 
            Pl_inds(3) = iP;
        elseif strcmp(Planets(iP).name,'Enceladus') 
            Pl_inds(4) = iP;
        end
    end
    plot_schematic(Planets(Pl_inds([1 2])),XLim);
    figure(367);set(gcf,'Name','Icy Planets Schematic');clf;
    set(gca,'XLim',[0 252],'Box','on','FontSize',20);
    plot_schematic(Planets(Pl_inds([3 4])),XLim);
%      for ij = 3:4         
%             R_m = Planets(Planet_inds(ij)).R_m;
%             T_oC = Planets(Planet_inds(ij)).T_oC;
%             P_MPa = Planets(Planet_inds(ij)).P_MPa;
%             depth_km = Planets(Planet_inds(ij)).depth_km;
%             z_cracking_km = Planets(Planet_inds(ij)).z_cracking_1mm_m/1e3;
%             d_ocean_km = Planets(Planet_inds(ij)).d_ocean_km;
%             
%                YLim = [Planets(Planet_inds(ij)).T_oC(1)-10 900];
%             
%             h=subplot('Position',[left bottom width height]);
%             line(depth_km,P_MPa,'Color',Planets(Planet_inds(ij)).plot_color,'LineStyle','--');%,[LineStyle{n}],'LineWidth',LineWidth)
%             t= text(0,P_MPa(1),[Planets(Planet_inds(ij)).name],'FontSize',HeadingSize);
%             set(gca,'YAxisLocation','left','XLim',XLim,'XDir','reverse');
%             ylabel('Pressure (MPa)','FontSize',LabelSize);
%            
%             ax1 = gca;
% 
%             ax2 = axes('Position',get(ax1,'Position'),...
%                     'YAxisLocation','right',...
%                     'Color','none',...
%                     'XColor','k','YColor','k');
%             line(depth_km,T_oC,'Color',Planets(Planet_inds(ij)).plot_color)
%             ylabel('Temperature (^oC)','FontSize',LabelSize);
%             set(gca,'YLim',YLim,'XLim',XLim,'XDir','reverse');
%             if ij == lP
%                 xlabel('Depth (km)','FontSize',LabelSize);
%             end
%             
%             box on;
%             bottom = bottom+height;
%             
%             subplot('Position',[left bottom width SchematicHeight])
%             vline(d_ocean_km,'b');
%             vline(z_cracking_km,'r');
%             set(gca,'XLim',XLim,'XDir','reverse','YTick',[],'XTick',[]);
%             box on;
%              
%             bottom = bottom-SchematicHeight-vert_space-3/2*height;
%     end
end
 %=========================================================================
function plot_schematic(Planets,XLim)
        lP = length(Planets);
        left = 0.1;
        height = 1/5;
        bottom = 1-1/10/2-3*height/2;
        vert_space = 0.3/lP;
        width = 0.8;
        SchematicHeight = height/2;
        LabelSize = 28;
        HeadingSize = 28;
        ArrowHeadStyle = 'vback1';
        AxisFontSize = 14;
        
    for ij = 1:lP         
            R_m = Planets(ij).R_m;
            T_oC = Planets(ij).T_oC;
            P_MPa = Planets(ij).P_MPa;
            depth_km = Planets(ij).depth_km;
            
            d_ocean_km = Planets(ij).d_ocean_km;
            z_cracking_km = Planets(ij).z_cracking_1mm_m/1e3 + d_ocean_km;
            
               YLim = [Planets(ij).T_oC(1)-10 900];
            
            h=subplot('Position',[left bottom width height]);
            line(depth_km,P_MPa,'Color',Planets(ij).plot_color,'LineStyle','--');%,[LineStyle{n}],'LineWidth',LineWidth)
            t= text(0,P_MPa(1),[Planets(ij).name],'FontSize',HeadingSize);
            set(gca,'YAxisLocation','left','XLim',XLim,'XDir','reverse','FontSize',AxisFontSize);
            if ij == lP
                ylabel('Pressure (MPa)','FontSize',LabelSize,'FontWeight','bold');
            end
            ax1 = gca;

            ax2 = axes('Position',get(ax1,'Position'),...
                    'YAxisLocation','right',...
                    'Color','none',...
                    'XColor','k','YColor','k');
            line(depth_km,T_oC(:,:,1),'Color',Planets(ij).plot_color)
            set(gca,'YLim',YLim,'XLim',XLim,'XDir','reverse','FontSize',AxisFontSize);
            if ij == lP
                xlabel('Depth (km)','FontSize',LabelSize,'FontWeight','bold');
                ylabel('Temperature (^oC)','FontSize',LabelSize,'FontWeight','bold');
                text(depth_km(end),P_MPa(end),'P','FontSize',20);
                text(depth_km(10),T_oC(1),'T','FontSize',20);
            end
            
            box on;
            bottom = bottom+height;
            
            subplot('Position',[left bottom width SchematicHeight])
            vline(d_ocean_km,'b');
            vline(z_cracking_km,'r');
            set(gca,'XLim',XLim,'XDir','reverse','YTick',[],'XTick',[],'FontSize',AxisFontSize);
            
            if ij==1
                text(0,0.2,'\leftarrowd\rightarrow');
                text(d_ocean_km(1),0.2,'\leftarrowz\rightarrow');
                
            end
            
            box on;
             
            bottom = bottom-SchematicHeight-vert_space-3/2*height;
    end
    end
%===============================================================================
function Planets = plot_tz(Planets)
global PLOT_ZT_1MM PLOT_ZT_10MM PLOT_H2 PLOT_ERRZ_VS_R
    lP = length(Planets);
            
if PLOT_ZT_1MM
	figure(370);
	for ij = 1:lP % plot cracking depth
		if Planets(ij).z_cracking_1mm_m
            [iy,ix] = size(Planets(ij).z_cracking_1mm_m);
            for jk = 1:iy
                z_km = 1e-3*Planets(ij).z_cracking_1mm_m(iy,:);
                if z_km>0
                    z_km = z_km(z_km>0);
                    t_Gyr = 1e-9*Planets(ij).t_yr(z_km>0);
                    plot(t_Gyr,z_km,...
                        'Color',Planets(ij).plot_color);%,[LineStyle{n}],'LineWidth',LineWidth)
                    t= text(t_Gyr(end),z_km(end),Planets(ij).name,...
                        'FontSize',18);
                end
                end
        end
    end
    if PLOT_H2
        figure(372); hold on
        for ij = 1:lP % plot dz/dt
            if Planets(ij).z_cracking_1mm_m
                z_m = Planets(ij).z_cracking_1mm_m;
                if z_m>0
                    z_m = z_m(z_m>0);
                    t_yr = Planets(ij).t_yr(z_m>0);
                    pp = spline(t_yr,z_m);
                    Planets(ij).dzdt_mm_yr = -1e3*pp1derv(pp,Planets(ij).t_yr); 
                    plot(Planets(ij).t_yr/1e6,Planets(ij).dzdt_mm_yr,'Color',Planets(ij).plot_color);%,[LineStyle{n}],'LineWidth',LineWidth)
                    t= text(t_yr(end)/1e9,Planets(ij).dzdt_mm_yr(end),Planets(ij).name,...
                        'FontSize',18);
                    ylabel('dz/dt (mm/yr)')
                end
            end
        end
    end
    if PLOT_ERRZ_VS_R
        figure(373)
        for ij = 1:lP
            R_km = 1e-3*Planets(ij).R_m;
            z_km = 1e-3*Planets(ij).z_cracking_1mm_m;
            if z_km>0
                z_km = z_km(z_km>0);
%                  errorbar(R_km,mean(z_km),z_km(end),z_km(1),...
%                    'Color',Planets(ij).plot_color)
                 errorbar(R_km,mean(z_km),z_km(end),z_km(1),...
                   'Color',Planets(ij).plot_color)
                text(R_km+10,mean(z_km),Planets(ij).name)
            end
        end
    end
end
if PLOT_ZT_10MM
	figure(371);
	for ij = 1:lP	
		if Planets(ij).z_cracking_10mm_m
            z_km = 1e-3*Planets(ij).z_cracking_10mm_m;
            if z_km>0
                z_km = z_km(z_km>0);
                t_Gyr = 1e-9*Planets(ij).t_yr(z_km>0);
        		plot(t_Gyr,z_km,'Color',...
                    Planets(ij).plot_color);%,[LineStyle{n}],'LineWidth',LineWidth)
    			t= text(t_Gyr(end),z_km(end),Planets(ij).name,'FontSize',18);
            end
		end
    end
    if PLOT_ERRZ_VS_R
        figure(374)
        for ij = 1:lP
            R_km = 1e-3*Planets(ij).R_m;
            z_km = 1e-3*Planets(ij).z_cracking_10mm_m;
            if z_km>0
                z_km = z_km(z_km>0);
%                  errorbar(R_km,mean(z_km),z_km(end),z_km(1),...
%                    'Color',Planets(ij).plot_color)
                 errorbar(R_km,mean(z_km),z_km(end),z_km(1),...
                   'Color',Planets(ij).plot_color)
                text(R_km+10,mean(z_km),Planets(ij).name)
            end
        end
    end
end
end
%==========================================================================
function text_rot(t_handle,x,y)
	global FontSize
    adj = (x(1)-x(end));
    opp = (y(1)-y(end));
    rot1 = atand(opp/adj); % inverse tangent, in degrees
    set(t_handle,'Rotation',-rot1);
end
%==========================================================================    
function [P_MPa,T_oC,rho_mc,z_m] = get_PT(d_ocean_m,R_planet_m,rho_av,H,To)
     global IS_EARTH lambda
     if IS_EARTH
      z_m = [0:d_ocean_m/10:d_ocean_m d_ocean_m+0.1:1e3:R_planet_m]; 
     else
               z_m = [0:1e3:R_planet_m]; 
     end
%    z_m = 0:10e3:(R_planet_m-d_ocean_m); %this works, but is slow
%z_m = 0:(R_planet_m-d_ocean_m)/200:(R_planet_m-d_ocean_m);

      
    mantle_inds = find(z_m>=d_ocean_m);
    ocean_inds = find(z_m<d_ocean_m);
    z_below_sf_m = z_m(mantle_inds)-d_ocean_m;
    
    T_oC(ocean_inds) = To;
    T_oC(mantle_inds) = To + z_below_sf_m./lambda * H./4/pi./(R_planet_m-d_ocean_m).^2;
 
    [P_MPa,rho_mc] = get_P(z_m,d_ocean_m,R_planet_m,rho_av);
end
%====================================================================
function [P_MPa,rho_mc] = get_P(z_m,d_ocean_m,R_planet_m,rho_av)
    % this is separate from get T because it is used by get_Pz_Cracking
    global IS_EARTH IS_MARS USE_P_EFFECTIVE P_PREM;
     mantle_inds = find(z_m>d_ocean_m);
    ocean_inds = find(z_m<=d_ocean_m);
    m2km = 1e-3;
    bar2MPa = 0.1;
	kbar2MPa = 1e2;
if IS_EARTH
   % [junk,junk,rho_mc,junk,junk,junk,P,junk] = PREM((R_planet_m-d_ocean_m-z_m)*1e-3);
    %P_PREM_MPa = kbar2MPa*P;
    P_MPa(ocean_inds) = 1e3*9.8*z_m(ocean_inds)/1e6;
	P_MPa(mantle_inds) = (1e3*9.8*d_ocean_m + 3500*10*(z_m(mantle_inds)-d_ocean_m))/1e6;
    rho_mc = 3500;
elseif IS_MARS
	P_MPa = 0 + 3500*6*z_m/1e6;
	rho_mc = 3500;
else
    if USE_P_EFFECTIVE
        [P_lith_bar,rho_mc] = d2P_2layer(R_planet_m*m2km,d_ocean_m*m2km,z_m*m2km,rho_av);
        P_hydr_bar = d2P_2layer(R_planet_m*m2km,(d_ocean_m+d_ocean_m)*m2km,z_m*m2km,rho_av);
        P_MPa = bar2MPa*(P_lith_bar-P_hydr_bar);
    else
        [P_bar,rho_mc] = d2P_2layer(R_planet_m*m2km,d_ocean_m*m2km,z_m*m2km ,rho_av);
        P_MPa = P_bar*bar2MPa;
    end
end
% OLD VERSION
% function [P,rho_mc] = get_P(z_m,d_ocean_m,R_planet_m,rho_av)
%     global IS_EARTH IS_MARS USE_P_EFFECTIVE P_PREM;
%     bar2MPa = 0.1;
% 	kbar2MPa = 1e2;
% if IS_EARTH
%     [junk,junk,rho_mc,junk,junk,junk,P,junk] = PREM((R_planet_m-d_ocean_m-z_m)*1e-3);
%     P_PREM = kbar2MPa*P;
%     inds = find(z_m>d_ocean_m)
%     P(1:length(z_m)-inds(1)-1) = 1e3*9.8*z_m/1e6;
% 	P(inds) = (1e3*9.8*d_ocean_m + 3500*10*z_m(inds))/1e6;
%     rho_mc = 3500;
% elseif IS_MARS
% 	P = 0 + 3500*6*z_m/1e6;
% 	rho_mc = 3500;
% else
%     if USE_P_EFFECTIVE
%         [P_lith,rho_mc] = d2P_2layer(R_planet_m,d_ocean_m,z_m+d_ocean_m,rho_av);
%         P_hydr = d2P_2layer(R_planet_m,d_ocean_m+d_ocean_m,z_m+d_ocean_m,rho_av);
%         P = bar2MPa*(P_lith-P_hydr);
%     else
%         [P,rho_mc] = d2P_2layer(R_planet_m,d_ocean_m,z_m + d_ocean_m,rho_av);
%         P = P*bar2MPa;
%     end
% end
end
%============================================================
 function deMartinCrackingDepthVGrainSize
global CrackingFront Tdot_str PLOT_H2 PLOT_ZT_10MM PLOT_ZT_1MM
% load cracking data for desired cooling rate and grain sizes as calculated
% using Calc_z_CrackingDeMartin2004.m

switch Tdot_str
        case '1oCyr'
        load('PcT10mm_1oCyr.ext');Pc10mm = PcT10mm_1oCyr(:,1);T10mm = PcT10mm_1oCyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
    	load('PcT1mm_1oCyr.ext'); Pc1mm = PcT1mm_1oCyr(:,1);T1mm = PcT1mm_1oCyr(:,2);% load Pc and T for Kic = 0.6 MPa m1/2
        load('PcTp1mm_1oCyr.ext');Pcp1mm = PcTp1mm_1oCyr(:,1);Tp1mm = PcTp1mm_1oCyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
        case '1oCMyr'
        load('PcT10mm_1oCMyr.ext');Pc10mm = PcT10mm_1oCMyr(:,1);T10mm = PcT10mm_1oCMyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
        load('PcT1mm_1oCMyr.ext'); Pc1mm = PcT1mm_1oCMyr(:,1);T1mm = PcT1mm_1oCMyr(:,2);% load Pc and T for Kic = 0.6 MPa m1/2
        load('PcTp1mm_1oCMyr.ext');Pcp1mm = PcTp1mm_1oCMyr(:,1);Tp1mm = PcTp1mm_1oCMyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
      case '1oCGyr'
        load('PcT10mm_1oCGyr.ext');Pc10mm = PcT10mm_1oCGyr(:,1);T10mm = PcT10mm_1oCGyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
        load('PcT1mm_1oCGyr.ext'); Pc1mm = PcT1mm_1oCGyr(:,1);T1mm = PcT1mm_1oCGyr(:,2);% load Pc and T for Kic = 0.6 MPa m1/2
        load('PcTp1mm_1oCGyr.ext');Pcp1mm = PcTp1mm_1oCGyr(:,1);Tp1mm = PcTp1mm_1oCGyr(:,2); % load Pc and T for Kic = 0.6 MPa m1/2
end

	CrackingFront = struct('T_1mm_oC',T1mm,'Pc_1mm_MPa',Pc1mm/1e6,...
		'T_10mm_oC',T10mm,'Pc_10mm_MPa',Pc10mm/1e6,...
		'T_p1mm_oC',Tp1mm,'Pc_p1mm_MPa',Pcp1mm/1e6);

if ~PLOT_H2
	figure(369);clf;hold on
    if PLOT_ZT_10MM
        plot(T10mm,Pc10mm/1e6,'k','LineWidth',3);
    end
    if PLOT_ZT_1MM
        plot(T1mm,Pc1mm/1e6,'k','LineWidth',3);
    end
	%plot(Tp1mm,Pcp1mm/1e6,'k','LineWidth',3);

	t=text(100,250,'Thermal Cracking Front: 1 mm','FontSize',16,'Rotation',14);
	%t= gtext('10 mm','FontSize',14,'Rotation',16);
	%t= gtext('0.1 mm','FontSize',14,'Rotation',16);
	text(200,50,'Cracking','FontSize',18,'FontWeight','bold');
	text(300,275,'No Cracking','FontSize',18,'FontWeight','bold');
end
 end
%==========================================================================
 function [P_cracking_MPa,T_cracking_oC,z_cracking_m] = get_Pz_cracking(P_front_MPa,T_front_oC,P_planet_MPa,T_planet,max_P_MPa,R_planet_m,d_ocean_m,max_z_m,rho_av)
        try
            pp_cracking_PT = spline(P_front_MPa,T_front_oC);
            pp_ocean_PT = spline(P_planet_MPa,T_planet);
			P_cracking_MPa = fzero(@(P_in) ppval(pp_cracking_PT,P_in)-ppval(pp_ocean_PT,P_in),[0 max_P_MPa]);
            T_cracking_oC = ppval(pp_cracking_PT,P_cracking_MPa);
            z_cracking_m = fzero(@(z) P_cracking_MPa-get_P(z,d_ocean_m,R_planet_m,rho_av),[d_ocean_m max_z_m]);
            z_cracking_m = z_cracking_m-d_ocean_m;
		catch
			P_cracking_MPa = 0;
            T_cracking_oC = 0;
            if R_planet_m < 500e3;
%     			z_cracking_m = R_planet_m-d_ocean_m;
                z_cracking_m = 0;
            else
                z_cracking_m = 0;
            end
     end
 end