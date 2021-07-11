%% Close and clear everything
close all;
clear all;
clc;
%% Hohmann Transfer code
MU=398600;
% Define keplerian elements of parking and target orbits
pA=9378;
pE=0;
pI=0;
pAsc=0;
pArgP=0;
pTA=0;
tA=35786;
tE=0;
tI=0;
tAsc=0;
tArgP=0;
tTA=0;

% Convert to Cartesian elements
[pR,pV]=coe2rv([pA pE pI pAsc pArgP pTA], MU);
[tR,tV]=coe2rv([tA tE tI tAsc tArgP tTA], MU);

pOrbitPeriod=floor(2*pi*sqrt(pA^3/MU));
targetOrbitPeriod=floor(2*pi*sqrt(tA^3/MU));
transferObritPeriod=floor(2*pi*sqrt((pA+tA)^3/(8*MU)));


eParkingOrbit=-MU/(2*pA);
eTransferObrit=-MU/(pA+tA);
eTargetOrbit=-MU/(2*tA);

vMagParkingOrbit=sqrt(2*(eParkingOrbit+(MU/pA)));
vMagTargetOrbit=sqrt(2*(eTargetOrbit+(MU/tA)));
vMagTransferPeriapsis=sqrt(2*(eTransferObrit+(MU/pA)));
vMagTransferApoapsis=sqrt(2*(eTransferObrit+(MU/tA)));

vTransferPeriapsis=[0;vMagTransferPeriapsis;0];
vTransferApoapsis=[0;vMagTransferApoapsis;0];

xinit=[pR' pV'];
% Integrate the nonlinear EoMs to get the points along the trajectory
% tolerance = 1e-6;
% options = odeset('RelTol',tolerance,'AbsTol',tolerance); 
% % tspan=[0 pOrbitPeriod-1];
% tspan=linspace(0,pOrbitPeriod-1,pOrbitPeriod);
% [pt,px] = ode45(@(t,x)two_body_eom(t,x,MU), tspan, xinit, options);
[pt,px] = rk4(1,0,pOrbitPeriod-1,xinit,MU);

xinit=[tR' tV'];
% Integrate the nonlinear EoMs to get the points along the trajectory
[tt,tx] = rk4(1,0,targetOrbitPeriod-1,xinit,MU);

xinit=[pR' vTransferPeriapsis'];
% Integrate the nonlinear EoMs to get the points along the trajectory
% [~,trx1] = rk4(1,0,targetOrbitPeriod/2,xinit,MU);
tspan=linspace(0,floor(transferObritPeriod/2),floor(transferObritPeriod/2)+1);
tolerance = 1e-8;
options = odeset('RelTol',tolerance,'AbsTol',tolerance); 
% Integrate the nonlinear EoMs to get the points along the trajectory
[~,trx1] = ode45(@(t,x)two_body_eom(t,x,MU), tspan, xinit, options);

xinit=[trx1(end,1:3) -vTransferApoapsis'];
tspan=linspace(floor(transferObritPeriod/2)+1,floor(transferObritPeriod),floor(transferObritPeriod/2)+1);
[~,trx2] = ode45(@(t,x)two_body_eom(t,x,MU), tspan, xinit, options);

figure 
plot3(px(:,1),px(:,2),px(:,3),'g','LineWidth',1.5);
hold on;
grid on;
plot3(tx(:,1),tx(:,2),tx(:,3),'r','LineWidth',1.5);
plot3(trx1(:,1),trx1(:,2),trx1(:,3),'m','LineWidth',1.5);
plot3(trx2(:,1),trx2(:,2),trx2(:,3),'y','LineWidth',1.5);
axis equal;

csvwrite('.\inner_orbit.csv',px);
csvwrite('.\outer_orbit.csv',tx);
csvwrite('.\transfer_orbit_to_outer.csv',trx1);
csvwrite('.\transfer_orbit_to_inner.csv',trx2);
