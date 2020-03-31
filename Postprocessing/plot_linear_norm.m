%plot linear norm

clear all
%close all

%gain
GT = [-0.1 0.1 -1];
%GT = [0.5];

%phys prop
visc = 0.5;
Ca = 6;

addpath(['~/Documents/MATLAB/droplet_simulations/drop_buoyancy/area_norm/lambda=' num2str(visc) '_Ca=' num2str(Ca)])

%time range
time = 0:0.1:10;

Area = zeros(numel(time),numel(GT));

figure
color = 'brgm';

%loop for different gains
for i = 1:numel(GT)
    
    gain = GT(i);
    
    filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta0.01_amp',num2str(gain),'_time0.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    Volume = axis_int_gauss(y',x');
    R = nthroot(3/4/pi*Volume,3);
    A0 = 4*pi*R^2;
    
    %loop for time
    for k = 1:numel(time)
        
        display(k)
        
        t = time(k);
    
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta0.01_amp',num2str(gain),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);
        
        Area(k,i) = surf_gauss_vect(y',x');
    
    end
    
    area_norm = (Area(:,i)-A0)/(Area(1,i)-A0);
    
    hold on
    plot(time,area_norm,color(i))
    grid on
    xlabel('t')
    ylabel('area variation norm')
    hold off
    
end