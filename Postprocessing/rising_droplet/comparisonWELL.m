%plot comparison between linear and non linear simulations at different
%time step

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

plotLIN = 0;
visc = 0.5;
Ca = 6;

%savinf option
save = 0;
dest = '~/Documents/phd_projects/viscous_drop_buoyancy/report14_compare_linear/figures/';

addpath(['~/Documents/MATLAB/droplet_simulations/drop_buoyancy/area_norm/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);
here = pwd;

folder = ['~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=' num2str(visc) '_Ca=' num2str(Ca)];
cd(folder)

amp = -1;
DELTA = 0.001;

load(['G=' num2str(amp) '_q=500_visc=' num2str(visc) '_Dt=0.001_loop=10000_DELTA='...
            num2str(DELTA) '_Ca=' num2str(Ca) '_RK2.mat'])

cd(here)
        
%time = 0:0.5:10;
time = 0:0.5:10;

%compute real area of the unperturbed sphere
a = risa(:,1)';
b = risb(:,1)';
Volume = axis_int_gauss(a,b);
R = nthroot(3/4/pi*Volume,3);
A0 = 4*pi*R^2;

AreaLinear = zeros(numel(time),1);
Vlinear = zeros(numel(time),1);
myarea = zeros(numel(time),1);
resizeV = zeros(numel(time),1);
resizeA= zeros(numel(time),1);
V_in = axis_int(risa(:,1)',risb(:,1)');
SAVEloop = zeros(numel(time),1);

%figure

for i=1:numel(time)
    
    disp(i)
    
    filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(D),'_amp-1_time',num2str(time(i)),'.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    loopNOW = round(time(i)/deltaT/checkpoint)+1;
    a = risa(1:q+1,loopNOW)';   b = risb(1:q+1,loopNOW)';
    
    xcm = center_mass(a,b);
    %linear data
    AreaLinear(i) = surf_gauss_vect(y',x');
    Vlinear(i) = axis_int_gauss(y',x');
    
    if plotLIN==1
    plot(x,y,-x,y)
    hold on
    axis equal
    drawnow
    grid on
    axis([-1.5 1.5 -1.5 1.5])
    end
    
    %there is an error of order delta^2 on the starting volume, compute
    %ratio between that and the analytical vualr
    resizeV(i) = (4*pi/3/Vlinear(i));
    resizeA(i) = resizeV(i)^(2/3);
    
    myarea(i) = Area(loopNOW);
    
    
end

plot(time,resizeV)
grid on

%compute nonlinear norm
AnormNONlinear = sqrt(Area-A0)/sqrt(Area(1)-A0);
AlinearResized = resizeA.*AreaLinear;
AnormLineaResized = sqrt(AlinearResized-4*pi)/sqrt(AlinearResized(1)-4*pi);

%load data drom expm calculation
load(['t0' num2str(visc*10) '' num2str(Ca)]);   t_expm = t;
load(['Anorm0' num2str(visc*10) '' num2str(Ca)]);   A_expm = Anorm;

figure
plot(0:checkpoint*deltaT:deltaT*loop,AnormNONlinear,'--k',t_expm,A_expm,time,AnormLineaResized,'LineWidth',2)
%plot(0:checkpoint*deltaT:deltaT*loop,AnormNONlinear,'--k','LineWidth',2)
hold on
hold off
%legend('Nonlinear','Linear expm','Linear direct A','Linear resized','Location','Best')
xlabel('t')
ylabel('||A||')
grid on
title(['\delta=' num2str(D)])
axis([0 time(end) 0.6 1.5])

if save==1
        
    name = 'manyCompare';
    print('-deps','-loose','-r100',[dest name sprintf('%04d',D*1000) '.eps'])
        
end

%plot linear volume evoulution
figure
plot(time,Vlinear-4*pi/3)
xlabel('t')
ylabel('V-V_{sphere}')
grid on