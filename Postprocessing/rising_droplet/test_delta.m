%test convesrion from elliptical perturbation to delta variation of the modes

close all
clear all

addpath('~/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED')

load('F_areaweight_mA1000.mat')

%phisical properties
visc = 5;
Ca = 6;
prol_obl = 0;

if prol_obl==1
    if visc==0.5 && Ca==6
        D = [0.031 0.032];
    end
    if visc==0.5 && Ca==5
        D = [0.045 0.046];
    end
    if visc==0.5 && Ca==4
        D = [0.07 0.08];
    end
    if visc==5 && Ca==6
        D = [0.048 0.049];
    end
    if visc==5 && Ca==5
        D = [0.069 0.07];
    end
    if visc==5 && Ca==4
        D = [0.103 0.104];
    end
    if visc==0.1 && Ca==6
        D = [0.042 0.043];
    end
    if visc==0.1 && Ca==5
        D = [0.06 0.061];
    end
    if visc==0.1 && Ca==7
        D = [0.031 0.032];
    end
elseif prol_obl==0
    if visc==0.5 && Ca==6
        D = -[0.083 0.084];
    end
    if visc==0.5 && Ca==5
        D = -[0.14 0.152];
    end
    if visc==0.5 && Ca==4
        D = -[0.19 0.24];
    end
    if visc==5 && Ca==6
        D = -[0.033 0.034];
    end
    if visc==5 && Ca==5
        D = -[0.052 0.053];
    end
    if visc==5 && Ca==4
        D = -[0.087 0.088];
    end
    if visc==0.1 && Ca==6
        D = -[0.16 0.3];
    end
    if visc==0.1 && Ca==5
        D = -[0.16 0.3];
    end
    if visc==0.1 && Ca==7
        D = -[0.25 0.255];
    end
end

%lambda=5 Ca=6
%D = -[0.0308 0.0319];
%D = [0.0478 0.048];
%D = -0.0319;

delta = zeros(1,numel(D));
deltaF = zeros(1,numel(D));

for i=1:numel(delta)
    
    display(i)
    [delta(i),deltaF(i),f] = find_delta(D(i),F);
    
end

% figure
% plot(D,delta,'ok-','Linewidth',2)
% xlabel('elliptic perturbation')
% ylabel('\delta')
% grid on
% title('delta naive')
% 
% figure
% plot(D,deltaF,'om-','Linewidth',2)
% xlabel('elliptic perturbation')
% ylabel('\delta')
% grid on
% title('delta area')

% figure
% plot([0 0],delta,'ok-','Linewidth',2)
% xlabel('elliptic perturbation')
% ylabel('\delta')
% grid on
% title('delta naive')
% 
% figure
% plot([0 0],deltaF,'om-','Linewidth',2)
% xlabel('elliptic perturbation')
% ylabel('\delta')
% grid on
% title('delta area')