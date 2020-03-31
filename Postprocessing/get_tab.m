%for unbounded poiseuille take converged cm velocity of different
%simulations

clear variables
close all

set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,'defaultlinelinewidth',3,'defaultpatchlinewidth',.7);

addpath('~/Documents/MATLAB/droplet_simulations/results/unbounded_poiseuille/')
%addpath('~/Documents/MATLAB/droplet_simulations/results/unbounded_poiseuille_2drops/2diff_size/')
addpath('~/Documents/MATLAB/droplet_simulations/server/')


%lambda = [0.1 0.5 1 5];
lambda = 1;
%Ca = [0.005 0.05 0.2 0.3 0.4 0.5];
Ca = 0.5;
alpha = 0.1:0.1:0.8;

%options
save = 1;
dest = '~/Documents/phd_projects/2droplets/report1_2unbounded/figures/';

%%%%%%%%%%%%%%%%%%%%%%%% PLOT 1 DROP %%%%%%%%%%%%%%%%%%%%%%

tab1 = zeros(numel(Ca),numel(alpha));
tab2 = zeros(numel(Ca),numel(alpha));

for i = 1:numel(lambda)
    for k = 1:numel(Ca)
        for j = 1:numel(alpha)
        
        filename = ['Poiseuille_Q=1_visc=' num2str(lambda(i)) '_Ca=' num2str(Ca(k)) '_R=1_alpha=' num2str(alpha(j)) '_RK2.mat'];
        load(filename)
        
        %give reaches convergence
        ind = find(risa(2,:)==0,1,'first');
        if isempty(ind)
            ind = PARAM.loop/PARAM.checkpoint;
        end
        xcm = zeros(ind-1,1);
        width = zeros(ind-1,1);
        
        for l = 1:ind-1
        %disp(l)
    
        xcm(l) = center_mass(risa(:,l)',risb(:,l)');
        width(l) = max(risb(:,l));
    
        end
        
        %v_max poiseuille
        v_max = 2*PARAM.Q/pi/PARAM.R^2;
        v_cm = diff(xcm(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;
        
        tab1(k,j) = v_cm(end);
        tab2(k,j) = width(end);
        
        end
        
    end
    
    figure(1)
    plot(alpha,tab1,'s--')
    xlabel('\alpha')
    ylabel('v_{cm}')
    title(['droplet velocity for \lambda=' num2str(lambda(i))])
    grid on
    
    figure(2)
    plot(alpha,tab2,'s--')
    xlabel('\alpha')
    ylabel('width')
    title(['droplet width for \lambda=' num2str(lambda(i))])
    grid on
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% PLOT 1 DROP %%%%%%%%%%%%%%%%%%%%%%

% Ca = 0.5;
% 
% tab1 = zeros(numel(Ca),numel(alpha));
% tab2 = zeros(numel(Ca),numel(alpha));
% 
% for i = 1:numel(lambda)
%     for k = 1:numel(Ca)
%         for j = 1:numel(alpha)
%         
%         filename = ['Poiseuille_Q=1_visc=' num2str(lambda(i)) '_Ca=' num2str(Ca(k)) '_R=1_alpha=' num2str(alpha(j)) '_RK2.mat'];
%         load(filename)
%         
%         %give reaches convergence
%         ind = find(risa(2,:)==0,1,'first');
%         xcm = zeros(ind-1,1);
%         width = zeros(ind-1,1);
%         
%         for l = 1:ind-1
%         %disp(l)
%     
%         xcm(l) = center_mass(risa(:,l)',risb(:,l)');
%         width(l) = max(risb(:,l));
%     
%         end
%         
%         %v_max poiseuille
%         v_max = 2*PARAM.Q/pi/PARAM.R^2;
%         v_cm = diff(xcm(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;
%         
%         tab1(k,j) = v_cm(end);
%         tab2(k,j) = width(end);
%         
%         end
%         
%     end
%     
%     figure(1)
%     hold on
%     plot(alpha,tab1,'ko-')
%     hold off
%     xlabel('\alpha')
%     ylabel('v_{cm}')
%     title(['droplet velocity for \lambda=' num2str(lambda(i))])
%     grid on
%     
%     figure(2)
%     hold on
%     plot(alpha,tab2,'ko-')
%     hold off
%     xlabel('\alpha')
%     ylabel('width')
%     title(['droplet width for \lambda=' num2str(lambda(i))])
%     grid on
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% PLOT 2 DROP %%%%%%%%%%%%%%%%%%%%%%

tab1 = zeros(numel(Ca),numel(alpha));
tab2 = zeros(numel(Ca),numel(alpha));
tab3 = zeros(numel(Ca),numel(alpha));

for i = 1:numel(lambda)
    for k = 1:numel(Ca)
        for j = 1:numel(alpha)
            
        disp(j)
        
%         if j<5
%             filename = ['2drops_el=' num2str(200*alpha(j)) '_Q=1_visc=1_Ca=0.05_Ca2=0.5_R=1_alpha1=' num2str(alpha(j)) '_alpha2=' num2str(alpha(j)) '_RK2.mat'];
%             load(filename)
%         else
%             filename = ['2drops_el=' num2str(400) '_Q=1_visc=1_Ca=0.05_Ca2=0.5_R=1_alpha1=' num2str(alpha(j)) '_alpha2=' num2str(alpha(j)) '_RK2.mat'];
%             load(filename)
%         end

        filename = ['2drops_Q=1_visc=1_Ca=' num2str(Ca) '_Ca2=' num2str(Ca) '_R=1_alpha1=' num2str(alpha(j)) '_alpha2=' num2str(alpha(j)) '_RK2.mat'];
        load(filename)
        
        %give reaches convergence
        ind = find(risa(2,:)==0,1,'first');
        if isempty(ind)
            ind = PARAM.loop/PARAM.checkpoint;
        end
        xcm = zeros(ind-1,1);
        width1 = zeros(ind-1,1);
        width2 = zeros(ind-1,1);
        dist = zeros(ind-1,1);
        
        for l = 1:ind-1
        %disp(l)
    
        xcm(l) = center_mass(risa(:,l)',risb(:,l)');
        width1(l) = max(risb(1:PARAM.q+1,l));
        width2(l) = max(risb(PARAM.q+2:end,l));
        dist = min(risa(1:PARAM.q+1,l)) - max(risa(PARAM.q+2:end,l));
    
        end
        
        %v_max poiseuille
        v_max = 2*PARAM.Q/pi/PARAM.R^2;
        v_cm = diff(xcm(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;
        
        tab1(k,j) = v_cm(end);
        tab2(k,j) = (width1(end) + width2(end))/2;
        tab3(k,j) = dist(end);
            
        end
        
    end
    
    %plot droplet velocity
    figure(1)
    hold on
    plot(alpha,tab1,'o-r')
    hold off
    xlabel('\alpha')
    ylabel('v_{cm}')
    title(['droplet velocity for \lambda=' num2str(lambda(i))])
    grid on
    legend(['1 drop Ca=' num2str(Ca)],'2 drops','Location','Best')
    
    if save==1
        
        name = ['vel_pair' num2str(Ca*100)];
        print('-deps','-loose','-r100',[dest name '.eps'])
        
    end
    
    %plot droplet width
    figure(2)
    hold on
    plot(alpha,tab2,'o-r')
    hold off
    xlabel('\alpha')
    ylabel('width')
    title(['droplet width for \lambda=' num2str(lambda(i))])
    grid on
    legend(['1 drop Ca=' num2str(Ca)],'2 drops','Location','Best')
    
    if save==1
        
        name = ['width_pair' num2str(Ca*100)];
        print('-deps','-loose','-r100',[dest name '.eps'])
        
    end
    
    %plot droplet distance
    figure(3)
    plot(alpha,tab3,'o-k')
    xlabel('\alpha')
    ylabel('dist')
    title(['droplet distance for \lambda=' num2str(lambda(i))])
    grid on
    
    if save==1
        
        name = ['dist' num2str(Ca*100)];
        print('-deps','-loose','-r100',[dest name '.eps'])
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%