%plot many droplets, given the data of a droplet evolution, plot droplet
%evolution for different timem (usually I plot the fisrt unstable)

close all
clear variables

Ca = 6;
visc = 0.5;

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
addpath(['~/Documents/MATLAB/droplet_simulations/drop_buoyancy/area_norm/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);
addpath(['~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);


%visc=0.5 Ca=6
%time = [0 3 3.3 6 9];
%time = [0 3 6.1 9 12];
%time = [0 1.15 3 6 9];
%time = [0 0.25 3 6 9];
%time = [0 3 6 9 12];
%time = [0 2.5 5 7.5 10];

%visc=5 Ca=6
%time = [0 7 14.4 21 28];
%time = [0 7 14 23.4 28];
%time = [0 7 14 25.6 32];
%time = [0 5.2 14 21 25];
%time = [0 7 14 21 30];

%visc=0.5 Ca=5
%time = [0 2.793 3 6 9];
%time = [0 2.903 3 6 9];
%time = [0 1.6 3 6 9];
%time = [0 0.27 3 6 9];

%visc=5 Ca=5
%time = [0.625 5 10 15 20];
%time = [0 10 15.691 20 30];
%time = [0 9.2 15 20 25];
%time = [1.567 10 15 20 23];

%color = 'mgbr';
color = 'bkkbrr';

%time = [0 2.5 5 7.5 10 14.5];
time = [0 2.5 5 7.5 10 14.5; 0 2.5 5 7.5 10 14.5; 0 2.5 5 7.5 10 14.5; 0 2.5 5 7.5 10 14.5; 0 2.5 5 7.5 10 14.5; 0 2.5 5 7.5 10 14.5];
amp = [-1 -1 -1 -1 -1 -1];
%D = [-0.0037 -0.0043 0.0061 -0.009];
DELTA = [-0.0140 -0.0135 -0.0140 0.0056 0.0055 0.0056];
%DELTA = [0.0061];

shift = 0;
vert = 0;

figure

for k = 1:numel(DELTA)
    
    this_color = color(k);
    %t = time(k);
    vert = 0;
    
    for i = 1:size(time,2)
        
        %vert = time(i);

        disp(i)
        
        if k==1 || k==4
        %plot LINEAR
        t = time(k,i);
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(DELTA(k)),...
            '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        hold on
        plot(x+shift,-y+vert,this_color,'Linewidth',2)
        plot(-x+shift,-y+vert,this_color,'Linewidth',2)
        axis equal
        hold off

        else
        
        %plot NONLINEAR
        t = time(k,i);
        load(['G=' num2str(amp(k)) '_q=200_visc=' num2str(visc) '_Dt=0.001_loop=15000_DELTA='...
            num2str(DELTA(k)) '_Ca=' num2str(Ca) '_RK2.mat'])
        
        ite = round(t/deltaT/checkpoint)+1;
        xcm = center_mass(risa(:,ite),risb(:,ite));
        hold on
        plot(-risb(:,ite)+shift,-risa(:,ite)+xcm+vert,this_color,'Linewidth',2)
        plot(risb(:,ite)+shift,-risa(:,ite)+xcm+vert,this_color,'Linewidth',2)
        axis equal
        hold off
        
        end

        %if i ~= 6
         %   vert = vert + time(k,i+1)-time(k,i);
        %end
        
        vert = vert + 2.5;

    end
    
    shift = shift + 2.5;
    %vert = vert - 2.5;

end

axis([-2 14 -2 14])
ylabel('t')
