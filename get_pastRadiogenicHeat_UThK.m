function H = get_pastRadiogenicHeat_UThK(M_planet_kg,f_Si,t_yr)
% Heat in W from radioactive decay of [U(238,235) Th(232) and K(40)]
% Co is a matrix containing the concentration (in ppm) of U, Th, and K
% Al26 is treated in Fish et al, 1960

%       [U     Th    K  ]
C_ppm = [0.012 0.040 840]*1e-6; % ppm, present concentrations of long-lived radionuclides in chondrites, Mason 1971
dEdC = [9.75 2.60 3.52e-4]*1e-5; % W kg-1, calc from Birth 1954 data
lambda = [1.551 0.495 5.543]*1e-10; % yr-1, Steiger and Jaeger 1977
for ij = 1:length(t_yr)
    H(ij) = M_planet_kg*f_Si*C_ppm*(dEdC.*exp(lambda.*t_yr(ij)))';
end  