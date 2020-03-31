%plot comparison between linear and non linear simulations at different
%time step

%clear all
close all

%load('~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/flow_field/ellipsoidal/q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.05_Ca=6_RK2.mat')

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,...
    'defaultpatchlinewidth',.7);

nFrames = 160;
step = 1;

scrsz = get(0,'ScreenSize');

%otptions
int = 1;    %plot interface
CM_vel = 0;    %plot center of mass velocity
save =1;    %save png
dest = '~/Documents/research_notes/APS/movie/dropOblateProlate/';

load('~/Documents/MATLAB/droplet_simulations/results/rising_dropletOLD/lambda=5_Ca=6/oblate/q=200_visc=5_Dt=0.005_loop=16000_DELTA=-0.14_Ca=6_RK2.mat')
risa1 = risa;   risb1 = risb;
load('~/Documents/MATLAB/droplet_simulations/results/rising_dropletOLD/lambda=5_Ca=6/prolate/q=200_visc=5_Dt=0.005_loop=16000_DELTA=0.06_Ca=6_RK2.mat')
risa2 = risa;   risb2 = risb;

%figure

for i = 3:nFrames
    
    ite = (i-1)*step+1;
    
    disp(i)
    
    if CM_vel==1
    
    a_before = risa(:,ite-2)';
    b_before = risb(:,ite-2)';
    a = risa(:,ite-1)';
    b = risb(:,ite-1)';

    %figure out center of mass velocity
    xcm_old = center_mass(a_before,b_before);
    xcm = center_mass(a,b);
    v_centermass = -(xcm-xcm_old)/deltaT/checkpoint;
    
    figure(2)
    hold on
    plot(i*checkpoint*deltaT,v_centermass,'o')
    hold off
    xlabel('t')
    ylabel('vcm')
    grid on
    
    end
    
    if int==1
        
    %load('~/Documents/MATLAB/droplet_simulations/results/rising_dropletOLD/lambda=5_Ca=6/oblate/q=200_visc=5_Dt=0.005_loop=16000_DELTA=-0.14_Ca=6_RK2.mat')
    
    a = risa1(:,ite)';
    b = risb1(:,ite)';
    
    xcm = center_mass(a,b);
    %a = a-xcm;
    
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
    
    %compute splines coordinates
    t = 0:0.05:0.9;
    ttt = repmat(t,1,numel(ax));
    axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
    bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
    cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
    dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
    ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
    byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
    cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
    dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
  
    %splines coordinates
    xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a(end)];
    yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b(end)];
    
    figure(1)
    plot(yyy+3.5,-xxx+xcm,'-b',-yyy+3.5,-xxx+xcm,'-b','LineWidth',2)
    axis equal
    axis off
    hold on
    axis([-1.5 5.5 -4.5 1.5])
    
%     if save==1
%         
%         name = 'oblate';
%         print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
% %         set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
% %         saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/interface' num2str(i) '.png'])
%     end
    
    %hold off
    
    %load('~/Documents/MATLAB/droplet_simulations/results/rising_dropletOLD/lambda=5_Ca=6/prolate/q=200_visc=5_Dt=0.005_loop=16000_DELTA=0.06_Ca=6_RK2.mat')
    
    a = risa2(:,ite)';
    b = risb2(:,ite)';
    
    xcm = center_mass(a,b);
    %a = a-xcm;
    
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
    
    %compute splines coordinates
    t = 0:0.05:0.9;
    ttt = repmat(t,1,numel(ax));
    axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
    bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
    cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
    dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
    ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
    byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
    cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
    dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
  
    %splines coordinates
    xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a(end)];
    yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b(end)];
    
    figure(1)
    hold on
    plot(yyy,-xxx+xcm,'-r',-yyy,-xxx+xcm,'-r','LineWidth',2)
    axis equal
    axis off
    axis([-1.5 5.5 -4.5 1.5])
    
    if save==1
        
        name = 'both';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
%         set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
%         saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/interface' num2str(i) '.png'])
    end
    
    hold off

    end
    
    
end

% Create AVI file.
%movie2avi(mov, 'myOblate.avi', 'compression', 'None');

