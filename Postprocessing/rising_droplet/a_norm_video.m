% compute area-variation norm taking into account the eventually different
% volume (from 4/3*pi)

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

dest = '~/Documents/research_notes/APS/movie/nonLinearEvo/';
save = 1;
plot2 = 1;
plotSeparate = 1;

%%%%%%%%%%%%%%%%%% LOAD NONLINEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/G=-1_q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.0056_Ca=6_RK2.mat')
risa1 = risa;   risb1 = risb;   Area1 = Area;
load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/G=-1_q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.0055_Ca=6_RK2.mat')
risa2 = risa;   risb2 = risb;   Area2 = Area;

range = 101;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD FOR LINEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ca = 6;
visc = 0.5;

addpath('/Users/Giacomo/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED')
addpath('/Users/Giacomo/Documents/MATLAB/stability/droplet_transient/subfunction')
%compute optimal initial shape
AAA = my_kojima(visc,Ca);     %kojima matrix
load(['finfo_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.mat'])
F = chol(Mwgt2); %weigth matrix
%I HAVE TO "CLEAN" fopt FROM delta!!!!!!!!!!
% [~,~,V] = svd(F*expm(0.2*AAA)/F);
% fopt1 = V(:,1);

T = time;
[~,j] = min(abs(2.95-T));
fopt = fopt(:,j);

%compute linear
DS = numel(1,range);
for i = 1:range
    
    disp(i)
    
    t = (i-1)*deltaT*checkpoint;
    
    %compute norm evoultion with matrix exponential
    f = expm(AAA*t)*fopt;
    DS(i) = norm(F*f)/norm(F*fopt);
    
end

figure(5)
plot((0:range-1)*deltaT*checkpoint,DS,'b')
xlabel('t')
ylabel('G_{\Delta S}')
axis([0 10 0 1.6])
grid on
drawnow

if save==1
        
        name = 'NonLinear0th';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
        
end

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

Area_norm1 = zeros(1,range);
Area_norm2 = zeros(1,range);
for i=1:range    
    
    disp(i)

    Area_norm1(i) = sqrt((Area1(i)-A01)/(Area1(1)-A01));
    Area_norm2(i) = sqrt((Area2(i)-A02)/(Area2(1)-A02));
    
    if plotSeparate==1
    
    figure(1)
    plot((0:range-1)*deltaT*checkpoint,DS,'b')
    hold on
    plot((0:i-1)*deltaT*checkpoint,Area_norm2(1:i),'k')
    plot((i-1)*deltaT*checkpoint,Area_norm2(i),'k*')
    grid on
    hold off
    xlabel('t')
    ylabel('G_{\Delta S}')
    axis([0 10 0 1.6])
    drawnow
    hold off
    
        if save==1

            name = 'NonLinear1st';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])

        end
    
        if plot2==1

        %plot area norm
        figure(4)
        plot((0:range-1)*deltaT*checkpoint,DS,'b')
        hold on
        plot((0:range-1)*deltaT*checkpoint,sqrt((Area2(1:range)-A02)/(Area2(1)-A02)),'k')
        plot((range-1)*deltaT*checkpoint,sqrt((Area2(range)-A02)/(Area2(1)-A02)),'*k')
        plot((0:i-1)*deltaT*checkpoint,Area_norm1(1:i),'r')
        plot((i-1)*deltaT*checkpoint,Area_norm1(i),'r*')
        grid on
        hold off
        xlabel('t')
        ylabel('G_{\Delta S}')
        axis([0 10 0 1.6])
        drawnow
        
        if save==1
        
        name = 'NonLinear2nd';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
        
        end

        end
    
    else
        
        figure(3)
        plot((0:range-1)*deltaT*checkpoint,DS,'b')
        hold on
        plot((0:i-1)*deltaT*checkpoint,sqrt((Area2(1:i)-A02)/(Area2(1)-A02)),'k')
        plot((i-1)*deltaT*checkpoint,sqrt((Area2(i)-A02)/(Area2(1)-A02)),'*k')
        plot((0:i-1)*deltaT*checkpoint,Area_norm1(1:i),'r')
        plot((i-1)*deltaT*checkpoint,Area_norm1(i),'r*')
        grid on
        hold off
        xlabel('t')
        ylabel('G_{\Delta S}')
        axis([0 10 0 1.6])
        drawnow
        
        if save==1

            name = 'NonLinearTogether';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])

        end
    
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
    plot(yyy1,-xxx1+xcm1,'r-',-yyy1,-xxx1+xcm1,'r-','LineWidth',3)
    axis equal
    axis([-1.3 1.3 -1.3 1.3])
    axis off
    
    if save==1
        
        name = 'dropUnstab';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
        
    end
    
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
    figure(3)
    plot(yyy2,-xxx2+xcm2,'k-',-yyy2,-xxx2+xcm2,'k-','LineWidth',3)
    axis equal
    axis([-1.3 1.3 -1.3 1.3])
    axis off
    
    if save==1
        
        name = 'dropStab';
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