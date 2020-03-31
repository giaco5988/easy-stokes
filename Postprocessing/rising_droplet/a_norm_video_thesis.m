% compute area-variation norm taking into account the eventually different
% volume (from 4/3*pi)

clear variables
close all

%set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

dest = '~/Documents/research_notes/presentationPHDdefense/movie/frames/';
save = 0;

%%%%%%%%%%%%%%%%%% LOAD NONLINEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathUP = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/forThesis/risingDroplet_transientGrowth/';
name1 = 'G=-0.1_q=200_visc=0.5_Dt=0.001_loop=10000_DELTA=0.0029_Ca=6_RK2.mat';
load([pathUP name1])
risa1 = risa;   risb1 = risb;   Area1 = Area;
risa1 = risa1(1:201,:); risb1 = risb1(1:201,:);
name2 = 'G=-1_q=200_visc=0.5_Dt=0.001_loop=10000_DELTA=0.0056_Ca=6_RK2.mat';
load([pathUP name2])
%load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/G=-1_q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.0055_Ca=6_RK2.mat')
risa2 = risa;   risb2 = risb;   Area2 = Area;

range = 101;
deltaT = PARAM.deltaT;
checkpoint = PARAM.checkpoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = risa1(:,1)';
b1 = risb1(:,1)';
a2 = risa2(:,1)';
b2 = risb2(:,1)';

Volume1 = axis_int_gauss(a1,b1);
Volume2 = axis_int_gauss(a2,b2);

R1 = nthroot(3/4/pi*Volume1,3);
R2 = nthroot(3/4/pi*Volume2,3);

%unperturbaed area
A01 = 4*pi*R1^2;
A02 = 4*pi*R2^2;

%Area = Area(1:range);

Area_norm1 = sqrt((Area1-A01)./(Area1(1)-A01));
Area_norm2 = sqrt((Area2-A02)./(Area2(1)-A02));
deltaS1 = Area1-A01;
deltaS2 = Area2-A02;
for i=1:range    
    
    disp(i)
    
    if real(Area_norm1(i))==0
        Area_norm1(i) = nan;
    end
    if Area_norm2(i)==0
        Area_norm2(i) = nan;
    end
    
    if i==1 && save==1
            width = 1201;
            height = 500;
            figure('pos',[50 100 width height])
    end
    
    figure(1)
    subplot(1,2,1)
    plot((0:i-1)*PARAM.deltaT*checkpoint,deltaS2(1:i),'r')
    hold on
    plot((i-1)*deltaT*checkpoint,deltaS2(i),'*r')
    plot((0:i-1)*deltaT*checkpoint,deltaS1(1:i),'k')
    plot((i-1)*deltaT*checkpoint,deltaS1(i),'k*')
    grid on
    hold off
    %xlabel('t')
    %ylabel('G_{\Delta S}')
    xyLabelTex('t','\Delta S')
    axis([0 10 0 0.025])
    %yticks(0:0.2:1.6)
    drawnow
        
%     if save==1
% 
%         name = 'NonLinearTogetherArea';
%         print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
% 
%     end
    
    figure(1)
    subplot(1,2,2)
    plot((0:i-1)*PARAM.deltaT*checkpoint,Area_norm2(1:i),'r')
    hold on
    plot((i-1)*deltaT*checkpoint,Area_norm2(i),'*r')
    plot((0:i-1)*deltaT*checkpoint,Area_norm1(1:i),'k')
    plot((i-1)*deltaT*checkpoint,Area_norm1(i),'k*')
    grid on
    hold off
    %xlabel('t')
    %ylabel('G_{\Delta S}')
    xyLabelTex('t','G_{\Delta S}')
    axis([0 10 0 1.6])
    yticks(0:0.2:1.6)
    drawnow
        
    if save==1

        name = 'NonLinearTogetherAll';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])

    end
    
    if i==1 && save==1
            width = 1201;
            height = 500;
            figure('pos',[50 100 width height])
    end
    
    %plot drop1 shape
    a1 = risa1(:,i);    b1 = risb1(:,i);
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a1', b1');
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
    xxx1 = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a1(end)];
    yyy1 = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b1(end)];
    xcm1 = center_mass(a1,b1);
    figure(2)
    plot(yyy1,-xxx1+xcm1,'k-',-yyy1,-xxx1+xcm1,'k-','LineWidth',3)
    %plot(b1,-a1+xcm1,'k-',-b1,-a1+xcm1,'k-','LineWidth',3)
    %axis equal
    %axis([-1.3 1.3 -1.3 1.3])
    %axis off
    
%     if save==1
%         
%         name = 'dropUnstab';
%         print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
%         
%     end
    
    %plot drop2 shape
    a2 = risa2(:,i);    b2 = risb2(:,i);
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a2', b2');
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
    xxx2 = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a2(end)];
    yyy2 = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b2(end)];
    xcm2 = center_mass(a2,b2);
    figure(2)
    hold on
    shift = 4;
    plot(yyy2+shift,-xxx2+xcm2,'r-',-yyy2+shift,-xxx2+xcm2,'r-','LineWidth',3)
    axis equal
    axis([-1.3 1.3+shift -1.3 1.3])
    axis off
    hold off
    drawnow
    
    if save==1
        
        %name = 'dropStab';
        name = 'twoDrops';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
        
    end

end

% figure(1)
% hold on
% plot((0:numel(Area)-1)*deltaT*checkpoint,Area_norm,'k')
% plot((0:5:numel(Area)-1)*deltaT*checkpoint,Area_norm(1:5:end),'sk')
% grid on
% xlabel('t')
% ylabel('Area variation norm')
% hold off