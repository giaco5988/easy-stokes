%compute nondimensional numbers for micromotor

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%real physical parameters (from Manjare et al. number 9 in TAB)
RealGamma = 72*1e-3;
RealRadius = 2*1e-6;
RealVisc =8.9*1e-4;
RealViscAir = 1.983*1e-5;
Mrate = 1.3*1e-14;
Qrate  = Mrate/1.2;
DeltaRho = 1e3-1.2;
g = 9.8;
rhoWater = 1e3;

%compute Capillary number
Ca = Qrate*RealVisc/RealRadius^2/RealGamma;
display(['Ca=' num2str(Ca)])

%compute Reynolds
Re = (rhoWater*Qrate)/(RealVisc*RealRadius);
display(['Re=' num2str(Re)])

%compute Bond
Bo = (DeltaRho*g*RealRadius^2)/RealGamma;
%Bo = (DeltaRho*g*RealRadius^4)/(RealGamma*RealVisc*Qrate);
%Bo = (DeltaRho*g*RealRadius^4)/(RealVisc*Qrate);
display(['Bo=' num2str(Bo)])

%compute viscosity ratio
lambda = RealViscAir/RealVisc;
display(['lambda=' num2str(lambda)])