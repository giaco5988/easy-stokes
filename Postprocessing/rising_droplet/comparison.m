%plot comparison between linear and non linear simulations at different
%time step

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%visc = 0.5 Ca=6
%time = [0 0.1 0.2  0.4 1.95 3.5 5.05 6.6 8.15 9.7];
%time = [0 0.2 1.44 3.92 7.64 10.12];
%time = 0:100:loop*deltaT;
%time = 0:end_time;
%time = [0 0.525 1.05 1.575 2.1 3.3 4.5 5.7 6.9 8.1 9.3];
%time = [0 1.325 2.65 3.975 5.3 6.625 7.95 9.275];
%time = [0];
%time = [0 0.825 1.65 2.475 3.3 4.125 5.775 6.6 7.425 8.25 9.075 9.9 10.725];
%time = [0 0.825 1.65 2.475 3.3 4.125 5.775 6.6 7.425];
%time = [0 0.825 1.65 2.475 3.3 4.125 5.775 6.6 7.425 8.25 9.075 9.9];
%time = [0 0.0625 0.125 0.1875 0.25 1.22 2.19 3.16 4.13 5.1 6.07 7.04 8.01 8.98]% 9.95];

%visc = 5 Ca=6
%time = [0 0.25 0.5 0.75 1 5.88 10.76 15.64 20.52 30.28 35.16];
%time = [0 3.6 7.2 10.8 14.4 18 21.6 25.2 28.8];

%visc=0.5 Ca=5
%time = [0 0.4 0.8 1.2 1.6 2.435 3.27 4.105 4.94 5.775 6.61 7.445 8.28 9.115];
%time = [0 0.4 0.8 1.2 1.6 2.435 3.27 4.105 4.94 5.775 6.61 7.445 8.28 9.115 9.95];

%visc=5 Ca=5
%time = [0 2.3 4.6 6.9 9.2 13.26 17.32 21.38 25.44 29.5 33.56 37.62];
%time = [0 2.3 4.6 6.9 9.2 13.26 17.32 21.38 25.44];

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
DELTA = -0.001;

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

a_lailai = zeros(numel(time),1);
v_lailai = zeros(numel(time),1);
elon_mine = zeros(numel(time),1);
elon_lailai = zeros(numel(time),1);
myarea = zeros(numel(time),1);
resizeV = zeros(numel(time),1);
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
    
%     elon_mine(i) = max(risa(1:nbrel(loop),loop))-min(risa(1:nbrel(loop),loop));
%     
%     elon_lailai(i) = max(y)-min(y);
    
    %xcm = center_mass(risa(1:nbrel(loop)+1,loop)',risb(1:nbrel(loop)+1,loop)');
    xcm = center_mass(risa(1:q+1,loopNOW)',risb(1:q+1,loopNOW)');
    a_lailai(i) = surf_gauss_vect(y',x');
    v_lailai(i) = axis_int_gauss(y',x');
    
    if plotLIN==1
    plot(x,y,-x,y)
    hold on
    axis equal
    drawnow
    grid on
    axis([-1.5 1.5 -1.5 1.5])
    end
    
    resizeV(i) = (V(loopNOW)*V_in/v_lailai(i))^(2/3);
    
    myarea(i) = Area(loopNOW);
    
    
end

%a_lailai(1) = Area(1);

b = a_lailai-A0;
bmod = a_lailai.*resizeV-A0;
myb = myarea-A0;

bNOmod = a_lailai-A0;

d = abs(elon_lailai-2);
myd = abs(elon_mine-2);

all_area = Area-A0;
a_norm = sqrt(all_area/all_area(1));

%load data drom expm calculation
load(['t0' num2str(visc*10) '' num2str(Ca)]);   t_expm = t;
load(['Anorm0' num2str(visc*10) '' num2str(Ca)]);   A_expm = Anorm;

figure
plot(0:checkpoint*deltaT:deltaT*loop,a_norm,'--k',t_expm,A_expm,time,sqrt(bNOmod/bNOmod(1)),time,sqrt(bmod/bmod(1)),'LineWidth',2)
hold on
hold off
legend('Nonlinear','Linear expm','Linear direct A','Linear resized','Location','Best')
xlabel('t')
ylabel('||A||')
grid on
%title(['\lambda=' num2str(visc) ' Ca=' num2str(Ca)])
title(['\delta=' num2str(D)])
axis([0 time(end) 0.6 1.5])

if save==1
        
    name = 'manyCompare';
    print('-deps','-loose','-r100',[dest name sprintf('%04d',D*1000) '.eps'])
        
end

figure
plot(time,v_lailai-4*pi/3)
xlabel('t')
ylabel('V')
grid on

% figure
% %plot((0:500:round(time(end)/deltaT))*deltaT,a_norm(1:500:round(time(end)/deltaT)+1),'-.sb',time,sqrt(bmod/bmod(1)),'k','LineWidth',2)
% plot(0:checkpoint*deltaT:deltaT*loop,a_norm,'-b',time,sqrt(bNOmod/bNOmod(1)),'k','LineWidth',2)
% hold on
% hold off
% legend('Nonlinear','Linear')
% xlabel('t')
% ylabel('Area norm')
% grid on
% title(['\lambda=' num2str(visc) ' Ca=' num2str(Ca)])
% axis([0 time(end) 1 1.5])
% 
% figure
% %plot((0:500:round(time(end)/deltaT))*deltaT,a_norm(1:500:round(time(end)/deltaT)+1),'-.sb',time,sqrt(bmod/bmod(1)),'k','LineWidth',2)
% plot(0:checkpoint*deltaT:deltaT*loop,a_norm,'-b',time,sqrt(bmod/bmod(1)),'k','LineWidth',2)
% hold on
% hold off
% legend('Nonlinear','Linear')
% xlabel('t')
% ylabel('Area norm resized')
% grid on
% title(['\lambda=' num2str(visc) ' Ca=' num2str(Ca)])
% axis([0 time(end) 1 1.5])