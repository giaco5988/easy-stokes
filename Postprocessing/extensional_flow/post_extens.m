%post processing extensional flow
clear variables
close all

%cd ~/Documents/MATLAB/droplet_simulations/results/val_extensional

addpath('~/Documents/MATLAB/droplet_simulations/results','~/Documents/MATLAB/droplet_simulations/results/extensional/viscosity/')
newton_dir = '~/Documents/MATLAB/droplet_simulations/results/extensional/newton_method/';
%newton_dir = '~/Documents/MATLAB/droplet_simulations/results/extensional/minimize/';

q = 100;
%allVisc = [0.1 1 10];
allVisc = 1;
%Ca = 0.08;

%allCa = 0.01:0.01:0.10;
%allq = 50:50:500;
%allq = 100;
%allD = zeros(1,numel(allCa));

%digitezed data
file_x1 = strcat('xxx_01.mat');
file_y1 = strcat('yyy_01.mat');
file_x2 = strcat('xxx_1.mat');
file_y2 = strcat('yyy_1.mat');
file_x3 = strcat('xxx_10.mat');
file_y3 = strcat('yyy_10.mat');

figure(4)


for k = 1:numel(allVisc)
    
    visc = allVisc(k);
    
    if visc==1
        allCa = 0.06:0.01:0.11;
    elseif visc==10
        allCa = 0.01:0.01:0.10;
    elseif visc==0.1
        allCa = 0.06:0.01:0.16;
    end
    
    allD = zeros(1,numel(allCa));
    %allDnewton = zeros(1,numel(allCa));


for i = 1:numel(allCa)
    
    disp(i)
    
    Ca = allCa(i);
    %q = allq(i);
    
    if visc == 10
    
    filename = strcat('Extensional_q=',num2str(q),'_visc=',num2str(visc),'_Dt=',num2str(0.005/Ca),'_loop=',...
            num2str(20000),'_Ca=',num2str(Ca),'_VR=',num2str(0),'_RK2.mat');
    
    if visc==10 && Ca==0.1
        filename = strcat('Extensional_q=',num2str(q),'_visc=',num2str(visc),'_Dt=',num2str(0.01/Ca),'_loop=',...
            num2str(10000),'_Ca=',num2str(Ca),'_VR=',num2str(0),'_RK2.mat');
    end
    
    elseif visc == 1
        
        filename = strcat('Extensional_q=',num2str(q),'_visc=',num2str(visc),'_Dt=',num2str(0.001/Ca),'_loop=',...
            num2str(10000),'_Ca=',num2str(Ca),'_VR=',num2str(0),'_RK2.mat');
        
    elseif visc == 0.1
        
        filename = strcat('Extensional_q=',num2str(q),'_visc=',num2str(visc),'_Dt=',num2str(0.02*Ca),'_loop=',...
            num2str(10000),'_Ca=',num2str(Ca),'_VR=',num2str(1),'_RK2.mat');
        
    end
    
    load(filename)
    
    L = (max(risa)-min(risa))/2;
    B = max(risb);
    D = (L-B)./(L+B);
    
%     figure(1)
%     hold all
%     plot(D(2:end)-D(1:end-1))
%     hold off
    
    ind = find(isnan(D),1,'first');

    %allD(i) = D(end);
    
    if visc==10 && Ca == 0.1 || visc == 0.1
        allD(i) = D(end);
    else
        allD(i) = D(ind-1);
    end
    
    figure(k)
    hold on
    plot(2:numel(D),D(2:end))
    grid on
    hold off
    title(['\lambda=' num2str(visc)])
    xlabel('ite')
    ylabel('D')
    
end

    %load results from newton method calculations
    try
    filename = [newton_dir 'RisingDroplet_n=50_visc=' num2str(visc) '.mat'];
    load(filename)
    allDnewton = manyD;
    allCanewton = manyCa;
    catch
    warning('no data')
    end

% file_x = strcat('xxx_',num2str(visc),'.mat');
% file_y = strcat('yyy_',num2str(visc),'.mat');
% 
% if visc == 0.1
% 
% file_x = strcat('xxx_01.mat');
% file_y = strcat('yyy_01.mat');
% 
% end

%load(file_x)
%load(file_y)

color = ['og-' 'or-' 'mo-'];

figure(4)
hold on
plot(allCa,allD,color(1+(k-1)*3:3*k),'LineWidth',2)
xlabel('Ca')
ylabel('D')
%legend('Leal data digitized','my DNS','Location','Best')
%title(strcat('Validation for \lambda=',num2str(visc)))
grid on
plot(allCanewton,allDnewton,'ok','LineWidth',2)
hold off

% figure
% plot(allq,allD,'o-')
% grid on
% xlabel('q')
% ylabel('D')
% title(strcat('Validation for  \lambda=',num2str(visc)))

end

figure(4)
hold on

load(file_x1);  load(file_y1);
plot(xxx,yyy,'--b','LineWidth',2)

if numel(allVisc)==3
%legend('my DNS \lambda=0.1','my DNS \lambda=1','my DNS \lambda=10','Stone at al. 1988','Location','Best')
legend('my DNS \lambda=0.1','my Newton \lambda=0.1','my DNS \lambda=1','my Newton \lambda=1','my DNS \lambda=10','my Newton \lambda=10','Stone at al. 1988','Location','Best')
elseif numel(allVisc)==2
legend('my DNS \lambda=1','my Newton \lambda=1','my DNS \lambda=10','my Newton \lambda=10','Stone at al. 1988','Location','Best')
end

load(file_x2);  load(file_y2);
plot(xxx,yyy,'--b','LineWidth',2)

load(file_x3);  load(file_y3);
plot(xxx,yyy,'--b','LineWidth',2)

xlabel('Ca')
ylabel('D')
%legend('Leal data digitized','my DNS','Location','Best')
%title(strcat('Validation for \lambda=',num2str(visc)))
hold off