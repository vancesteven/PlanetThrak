function f = get_SiMassFraction(rho_ave,rho_ice,rho_Si)
% get silicate mass fraction of an icy body as per Schubert et al
% 1986, page 230
f = (1-rho_ice./rho_ave)./(1-rho_ice./rho_Si);