function plotMicrofractureStress
% this is an older file written prior to reproducing the model's
% details.  It's a sound idea to use a linearized model, but it's
% reassuring that the results could be reproduced.

% Reproduction of the critical results from DeMartin Et Al 2004.
% Linearization of P and T at which microfracturing occurs in olivine,
% assuming this boundary is defined by critical stress intensity, K_IC =
% 0.6 MPa m^(1/2), as measured by the authors (+- 0.3) for olivine, and modeled
% using Friedrich and Wong (1986) and Evans and Clark (1980).

figure(3333);clf
delT = 0:100:1200;

b = 500;
m = (0-160)/(b-1200);
P_100micron = m*(delT - b);

b = 300;
m = (50-300)/(b-1100);
P_1mm = m*(delT - b);

plot(delT,P_100micron,delT,P_1mm,'-.');
axis([0 1200 0 300]);
set(gca,'YDir','reverse');
title('Thermal Cracking');
ylabel('Pressure (MPa)');
xlabel('\DeltaT (K)');
gtext('100 \mum');
gtext('1 mm');