%compute delta prime using the F for the area variation norm

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variable
visc = 0.5;
Ca = 6;

addpath('~/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED');
addpath(['~/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED/lambda='...
    num2str(visc) '_Ca=' num2str(Ca)]);

%delta values before rescaling, choose if delta stable or unstable
filedelta = ['delta_lambda' num2str(visc) '_Ca' num2str(Ca) '.mat'];
load(filedelta);
delta_low = delta(1,:);
delta_high = delta(2,:);
delta = sum(delta)/2;

%initialization
delta_star1 = zeros(1,numel(delta_low));
delta_star2 = zeros(1,numel(delta_high));
delta_star = zeros(1,numel(delta));

% load data
filename1 = ['finfo_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.mat']; %CAREFUL THERE IS A TIME HERE!!!
load(filename1);
filename3 = ['Gt_area_t_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.dat'];
A = importdata(filename3);
t = A(:,1)';
filename2 = ['time_lambda' num2str(visc) '_Ca' num2str(Ca) '.mat'];
load(filename2);
%filename4 = ['prolate_lambda' num2str(visc) '_Ca' num2str(Ca) '.mat'];
%load(filename4);
%filename5 = ['oblate_lambda' num2str(visc) '_Ca' num2str(Ca) '.mat'];
%load(filename5);

%load F
load('F_areaweight_mA1000.mat');
%M = Mwgt2;
%F = chol(M);

for iii = 1:numel(delta)

    [~,j] = min(abs(time(iii)-t));

    %initial shape
    f = fopt(:,j);

    %correction factor CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %new_norm = norm((F*f)'*(F*f));
    new_norm = norm((F*f));
    %new_norm = sqrt((F*f)'*(F*f));

    %delta prime
    delta_star(iii) = delta(iii)*new_norm;
    delta_star1(iii) = delta_low(iii)*new_norm;
    delta_star2(iii) = delta_high(iii)*new_norm;
    
    figure(1)
    hold on
    plot([time(iii) time(iii)], [delta_low(iii) delta_high(iii)],'o-k')
    hold off
    
    figure(2)
    hold on
    plot([time(iii) time(iii)], [delta_star1(iii) delta_star2(iii)],'o-k')
    hold off
    
    figure(4)
    hold on
    plot([time(iii) time(iii)], [delta_star1(iii) delta_star2(iii)],'o-k')
    hold off
    
end

figure(1)
hold on
plot(time,delta,'o-')
xlabel('T')
grid on
ylabel('\delta')
hold off
title(['\lambda=' num2str(visc) ' and Ca=' num2str(Ca)])

% filename4 = ['~/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED/lambda=' num2str(visc) '_Ca=' num2str(Ca)];
% cd(filename4)
% savefig('old_delta')

figure(2)
hold on
plot(time,delta_star,'o-r')
xlabel('T')
ylabel('\delta^*')
grid on
plot([time(1) time(7)],sum(delta_pro)/2*ones(2,1),'--b')
plot(ones(2,1)*time(7),delta_pro,'ok-')
plot([time(1) time(7)],sum(delta_obl)/2*ones(2,1),'--g')
plot(ones(2,1)*time(7),delta_obl,'ok-')
hold off
title(['\lambda=' num2str(visc) ' and Ca=' num2str(Ca)])

%savefig('new_delta')

% figure
% plot(A(:,1),A(:,2))
% xlabel('T')
% ylabel('Gt')
% grid on
% title(['transient growth \lambda=' num2str(visc) ' and Ca=' num2str(Ca)])

[~,i] = min(abs(time(end)-A(:,1)));

figure(4)
hold on
[ax,line1,line2] = plotyy(time,delta_star,A(1:i,1),A(1:i,2));
plot([0 time(7)],min(sum(delta_pro),sum(delta_obl))/2*ones(2,1),'--k')
plot(ones(2,1)*time(7),[min(delta_pro(1),delta_obl(1)) min(delta_pro(2),delta_obl(2))],'ok-')
hold off
grid on
xlabel('T')
ylabel(ax(1),'\delta_{crit}')
ylabel(ax(2),'G_{\Delta S}')
title(['\lambda=' num2str(visc) ' and Ca=' num2str(Ca)])

set(ax(1),'xlim')
set(line1,'Color','k')
set(line1,'Marker','o')
set(line2,'Color','k')
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','k')
set(ax(1),'xlim',[0 26])
set(ax(1),'ylim',[0.048 0.056])
set(ax(2),'ylim',[1 1.2])
%set(ax(1),'ylim',[0.049 0.06])
set(ax(2),'xlim',[0 26])
%text(3.5,0.05,'\delta_{crit} prolate')
set(ax(1),'YTick',0.048:0.002:0.056)
set(ax(2),'YTick',1:0.05:1.2)



%savefig('transient_growth')