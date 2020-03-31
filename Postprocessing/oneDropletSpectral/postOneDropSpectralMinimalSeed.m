%post processing of test time stepping spectral

close all
clear variables

%data
Ca = 0.1;
visc = 1;
n = 25;
Tend = 10;
BC = 1;         % 1 is extensional flow, 2 is rising droplet
legendre = 1;
Thorizon = 20;
EnergyPerturb = 0.1;
res = 0;        %results or server

%options
plotShape = 1;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
end

%filename
name = ['minimalSeedSpectralXYmodes_n=' num2str(n) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_T=' num2str(Thorizon) '_A0=' num2str(EnergyPerturb) '.mat'];

%upload data
upload = [dir name];
load(upload)

%PARAM.legendre = 0;

%initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

T = T{end};
Y = Y{end};

%loops
AreaPost = zeros(1,numel(T));
for i = 1:numel(T)
    
    xMode = Y(i,1:2:end-1)';
    yMode = Y(i,2:2:end)';
    [x,y] = fromModesToGrid(xMode,yMode,PARAM);
    
    AreaPost(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM);
    
    if plotShape==1
        
        figure(1)
        plot(x,y,'k')
        hold on
        plot(x,-y,'k')
        axis equal
        axis([-2 2 -2 2])
        grid on
        xlabel('z')
        ylabel('r')
        title('drop shape')
        hold off
        drawnow
        
    end
    
end

%plot initial and final shape
figure
plot([xBase; flip(xBase)],[yBase; -flip(yBase)],'k')
hold on
plot([x; flip(x)],[y; -flip(y)],'r')
axis equal
axis ([-2 2 -2 2])
xlabel('z')
ylabel('r')
%title('droplet shape')
legend('Initial shape','Final shape','Location','Best')
grid on
drawnow

%plot objective function
figure
plot(-OutFun0,'x-')
hold on
%plot(areaBaseState+PARAM.A0perturb*ones(1,PARAM.iterMinimalSeed))
xlabel('ite')
ylabel('F')
title('Optimization')
grid on

%plot surface area
AreaNorm = (AreaPost-areaBaseState)/(AreaPost(1)-areaBaseState);
figure
plot(T,AreaNorm)
grid on
xlabel('T')
ylabel('||\Delta S||')
title('Surface area norm')

%print simulation time
display(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











