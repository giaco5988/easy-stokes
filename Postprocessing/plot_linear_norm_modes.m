%plot linear norm

clear all
%close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

addpath('/Users/Giacomo/Documents/MATLAB/stability/droplet_transient/subfunction')
addpath('/Users/Giacomo/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED')

%gain
%GT = [-0.1 0.1 -1];
%GT = [0.5];

%time at which I maximize the modes
%t_max = [0.625 15.691 9.2 1.567];
%t_max = [0.25 1.25 3.3 6.1];
%t_max = [5.2 14.4 23.4 25.6];
t_max = 3.3;

%phys prop
visc = 0.5;
Ca = 6;
t_range = 0:0.2:10;

%kojima matrix
A = my_kojima(visc,Ca);

%load initial shapes f and F
load('F_areaweight_mA1000')
load(['finfo_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.mat'])   %CAREFUL THERE IS A TIME HERE!!!
AAA = importdata(['Gt_area_t_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.dat']);

color = 'mgbr';

t = AAA(:,1)';
Gt = AAA(:,2)';

figure
plot(t,Gt,'k')

%loop for different maximization
for i = 1:numel(t_max)

    [~,j] = min(abs(t_max(i)-t));

    %initial shape
    f = fopt(:,j);
    %load('f_oblate')
    %load('f_oblate_056')
    
    %time evolution
    A_norm = zeros(numel(t_range),1);
    
    for k = 1:numel(t_range)
        
        display(k)
        A_norm(k) = norm(F*expm(A*t_range(k))*f)/norm(F*f);
        %A_norm(k) = norm(F*expm(A*t_range(k))*inv(F));
        
    end
    
    hold on
    plot(t_range,A_norm,color(i))
    %draw circle on the envelope
    %plot(t(j),Gt(j),[color(i) 'o'])
    hold off

end

grid on
legend('envelope',['t_{max}=' num2str(t_max(1))],['t_{max}=' num2str(t_max(2))],['t_{max}=' num2str(t_max(3))],['t_{max}=' num2str(t_max(4))])
title(['\lambda=' num2str(visc) ' Ca=' num2str(Ca)])
xlabel(['$ t $'],'interpreter','latex')
ylabel(['$ G_{\Delta S} $'],'interpreter','latex')

filename = ['~/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED/lambda=' num2str(visc) '_Ca=' num2str(Ca)];
cd(filename)
%savefig('transient_growth')

% print('-depsc','-r300',['~/Documents/Phd_projects/Viscous_drop_buoyancy/Report8_allgraphs/lambda' num2str(visc) '_Ca' num2str(Ca)...
%         '/figures/transient_growth_56.eps'])