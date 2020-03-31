%compare drop shapes, linear and non linear

close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

addpath(['~/Documents/MATLAB/droplet_simulations/drop_buoyancy/area_norm/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);

%time when I compare
%visc=0.5 Ca=5
%time = [0.4 1.6 3.27 5.775 8.28 9.95];
%visc=0.5 Ca=6
%time = [0.825 3.3 5.775 6.6 9.9];
%visc=5 Ca=6
time = [3.6 14.4 21.6 25.2 28.8];

%shift nonlinear droplet
shift = 2.5;

%loop for different times
for i = 1:numel(time)
    
    display(i)
    
    filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(D),'_amp-1_time',num2str(time(i)),'.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    figure
    plot(x,-y,'r')
    hold on
    plot(-x,-y,'r')
    
    ite = round(time(i)/deltaT)+1;
    a = risa(:,ite)';
    b = risb(:,ite)';
    xcm = center_mass(a,b);
    
    plot(b+shift,-a+xcm,'b')
    plot(-b+shift,-a+xcm,'b')
    
    axis equal
    axis([-1.2 3.7 -1.5 1.5])
    title(['t=' num2str(time(i))])
    
    %print('-depsc','-r300',['~/Documents/Phd_projects/Viscous_drop_buoyancy/Report8_allgraphs/lambda' num2str(visc) '_Ca' num2str(Ca)...
     %   '/figures/drops_t' num2str(1000*time(i)) '.eps'])
    
end