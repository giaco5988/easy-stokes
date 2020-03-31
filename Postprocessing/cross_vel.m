%compute velocity of the center of mass and of the tip of the droplets

close all
clear all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical prpoerties
visc = 5;
Ca = 5;

%choose color of the graph
color = 'rkgmc';
figure

filename1 = ['~/Documents/MATLAB/droplet_simulations/results/delta_graph/lambda=' num2str(visc) '_Ca=' num2str(Ca) '/optimal/'];
cd(filename1)
%addpath(['~/Documents/MATLAB/droplet_simulations/results/delta_graph/lambda=' num2str(visc) '_Ca=' num2str(Ca) '/optimal/'])

%geometrical and file properties
%DELTA = [0.0085 0.0086 0.0087];
DELTA = 0.0086;
Gt = -1;
dt = 0.02;
end_loop = 4000;
maxk = 1400;

%conuter
w = 1;

%loop for delta
for i = 1:numel(DELTA)
    
    delta = DELTA(i);
    filename2 = ['G=' num2str(Gt) '_q=200_visc=' num2str(visc) '_Dt=' num2str(dt) '_loop=' num2str(end_loop)...
        '_DELTA=' num2str(delta) '_Ca=' num2str(Ca) '_RK2.mat'];
    load(filename2)
    
    xcm_old = center_mass(risa(:,1)',risb(:,1)');
    
    %figure
    hold all
    for k = 1:end_loop
        
        disp(k)
        
        a = risa(:,k)';
        b = risb(:,k)';
        
        xcm = center_mass(a,b);
        v_centermass = -(xcm-xcm_old)/deltaT;
        v_tip = -risy(1,k);
        
        if k==maxk||v_tip==0
            break
        end
        
        if abs((v_tip-v_centermass)/v_centermass) < 0.0001
            ite_cross(w) = k;
            w = w+1;
            %break
        end
        
        if k>1
            %hold on
            plot(k*deltaT,v_centermass,'b')
            plot(k*deltaT,v_tip,color(i))
            %hold off
        end
        
        xcm_old = xcm;
        
    end
    
    grid on
    xlabel('t')
    ylabel('v')
    %legend('v_{centermass}','v_{tip}','Location','Best')
    title(['\lambda=' num2str(visc) ' Ca=' num2str(Ca) ' \delta=' num2str(DELTA)])
    
%     print('-depsc','-r300',['~/Documents/Phd_projects/Viscous_drop_buoyancy/Report8_allgraphs/lambda' num2str(visc) '_Ca' num2str(Ca)...
%         '/figures/crossv_delta' num2str(delta) '.eps'])
    
end

%print('-depsc','-r300',['~/Documents/Phd_projects/Viscous_drop_buoyancy/Report8_allgraphs/lambda' num2str(visc) '_Ca' num2str(Ca)...
        %'/figures/crossv_moredelta' num2str(delta) '.eps'])