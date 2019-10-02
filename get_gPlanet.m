function g_planet = get_gPlanet(d,D_ocean_km,R_planet_km,rho_av_kgm3)
% d and D_ocean_km are in km
run V_consts;

%rho_rock = 3800;rho_water=1200; % kg/m3  -- rho_rock is chosen based on Anderson et al (1997), where they indicate from C/MR2 data that the rock density for a Europa with no core is the chosen value
rho_rock = rho_m2Layer(R_planet_km,D_ocean_km,rho_av_kgm3);
rho_water = 1000;
if find(d>=D_ocean_km)
    M_ocean(find(d>=D_ocean_km)) = 0;
    M_ocean(find(d<D_ocean_km)) = ((D_ocean_km-d(find(d<D_ocean_km)))*km2m).^3*rho_water;
    R_mantle(find(d<D_ocean_km)) = (R_planet_km - D_ocean_km)*km2m;
    R_mantle(find(d>=D_ocean_km)) = (R_planet_km-d(find(d>=D_ocean_km)))*km2m;
else M_ocean = ((D_ocean_km-d)*km2m).^3*rho_water;R_mantle = (R_planet_km - D_ocean_km)*km2m;
end
M = 4/3*pi*((R_mantle).^3*rho_rock+M_ocean);
g_planet = G*M./((R_planet_km-d)*km2m).^2;