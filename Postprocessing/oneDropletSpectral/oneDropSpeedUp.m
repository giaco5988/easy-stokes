%compare results of test time stepping spectral

close all
clear variables

%data
manyCa = 6*ones(4,1);
visc = 1*ones(4,1);
n = 50*ones(4,1);
manyTend = 20*ones(16,1);
manyODE = [1 1 1 1];                        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 2*ones(16,1);                         % 1 is extensional flow, 2 is rising droplet
volCorr = 0*ones(1,4);
MANYmaxDT = 1e-2*ones(4,1);
legendre = 1*ones(4,1);
manyD = 0.01*ones(4,1);
manyCPUs = 1:4;

%number of comparison
compare = numel(manyCPUs);

%option
plotShape = 0;
plotModes = 0;
plotRes = 1;
Tbreak = 200;
step = 1;

%directory
%dir = '~/Documents/MATLAB/droplet_simulations/server/';
dir = '~/Documents/MATLAB/droplet_simulations/results/';

%filename
%name = ['DropSpectral_ODE=' num2str(manyODE(1)) '_Legendre=' num2str(legendre(1)) '_BC=' num2str(BC(1)) '_Ca=' num2str(manyCa(1)) '_visc=' num2str(visc(1)) '_n=' num2str(n(1)) '_maxDT=' num2str(MANYmaxDT(1)) '_volCorr=' num2str(volCorr(1)) '_Tend=' num2str(manyTend(1)) '.mat'];
name = ['DropSpectral_ODE=' num2str(manyODE(1)) '_Legendre=' num2str(legendre(1)) '_BC=' num2str(BC(1)) '_Ca=' num2str(manyCa(1)) '_visc=' num2str(visc(1)) '_n=' num2str(n(1)) '_D=' num2str(manyD(1)) '_maxDT=' num2str(MANYmaxDT(1)) '_volCorr=' num2str(volCorr(1)) '_Tend=' num2str(manyTend(1))  '_CPUs=' num2str(manyCPUs(1)) '.mat'];

%upload data
upload = [dir name];
load(upload)

%initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%initialize
time = zeros(1,compare);

%loop for comparison
for l = 1:compare
    
%print to screen
display([num2str(l) ' of ' num2str(compare)])
    
%filename
name = ['DropSpectral_ODE=' num2str(manyODE(l)) '_Legendre=' num2str(legendre(l)) '_BC=' num2str(BC(l)) '_Ca=' num2str(manyCa(l)) '_visc=' num2str(visc(l)) '_n=' num2str(n(l)) '_D=' num2str(manyD(l)) '_maxDT=' num2str(MANYmaxDT(l)) '_volCorr=' num2str(volCorr(l)) '_Tend=' num2str(manyTend(l)) '_CPUs=' num2str(manyCPUs(l)) '.mat'];

%upload data
upload = [dir name];
load(upload)

time(l) = simulationTime;

end

%speed up
figure
speedUp = time(1)./time;
plot(speedUp,'o-')
hold on
plot(manyCPUs,manyCPUs)
grid on
xlabel('CPUs')
ylabel('Speed Up')
title('Speed up in Matlab')











