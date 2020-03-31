%micromotor evoultion till convergence (hopefully)

close all

%number of elements
n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;

%options
plotRT = 0;
compareData = 1;
PARAM.continuation = nan;

%give reaches convergence
ind = find(risa(2,:)==0,1,'first');
if isempty(ind)==1
    ind = PARAM.loop/PARAM.checkpoint;
end

%ind = PARAM.loop/PARAM.checkpoint;

%volume preallocation
V = zeros(ind,1);
Fdrop = zeros(ind,1);

%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i = 1:ind-1
    
    disp(i)
    
    %wall coordinates
    xWall = risa(1:n+m+j+1,i);  yWall = risb(1:n+m+j+1,i);
    
    %droplet coordinates
    xDrop = risa(n+m+j+2:n+m+j+q+2,i);   yDrop = risb(n+m+j+2:n+m+j+q+2,i);
    
    %plot drop shape
    if plotRT==1
        figure(1)
        plot(xWall,yWall,'k-',xWall,-yWall,'k-')
        hold on
        plot(xDrop,yDrop,'r-',xDrop,-yDrop,'r-')
        title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ' \alpha=' num2str(PARAM.alpha)])
        xlabel('x')
        ylabel('r')
        hold off
        grid on
        axis equal
        drawnow
        hold off
    end
    
    %drop volume
    V(i) = axis_int_gauss_vect(xDrop',yDrop');
    
    %force actinf on droplet
    [df_x,df_y,K,K1,K2] = df_surf_tens_spline_symm(xDrop',yDrop',1);
%     K1 = curv_spline(xDrop',yDrop');
%     K2 = N(2,:)./yDrop';
%     K2(1) = K1(1);
%     K2(end) = K1(end);
    Fdrop(i) = int_axis_spline_symmetric(xDrop',yDrop',df_x');
    
    %dl = sqrt(diff(xDrop).^2+diff(yDrop).^2);
    %Fdrop(i) = pi*sum(dl'.*(df_x(1:end-1).*yDrop(1:end-1)'+df_x(2:end).*yDrop(2:end)'));
    
end

figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

%wall coordinates
xWall = risa(1:n+m+j+1,ind-1);  yWall = risb(1:n+m+j+1,ind-1);
    
%droplet coordinates
xDrop = risa(n+m+j+2:n+m+j+q+2,ind-1);   yDrop = risb(n+m+j+2:n+m+j+q+2,ind-1);

%plot converged droplet shape
subplot(2,2,1)
plot(xWall,yWall,'k-',xWall,-yWall,'k-')
hold on
plot(xDrop,yDrop,'r-',xDrop,-yDrop,'r-')
title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ' \alpha=' num2str(PARAM.alpha)])
xlabel('x')
ylabel('r')
grid on
hold off
axis equal

%plot velocity residuals
subplot(2,2,2)
try
    semilogy((2:ind-1)*PARAM.checkpoint,residuals(2:ind-1),'LineWidth',2)
catch exception
    semilogy((2:ind-1)*PARAM.checkpoint,many_res(2:ind-1),'LineWidth',2)
end
grid on
xlabel('iteration')
ylabel('res')
title('velocity residuals')

%plot curvature in the final configuration
% subplot(2,2,3)
% plot(risa(n+m+j+2:end,ind-1),KKK(:,ind-2),'o-')
% title('curvature')
% xlabel('x')
% ylabel('k')
% grid on

%plot droplet velocity
subplot(2,2,3)
Umax = 2*PARAM.Q/pi/PARAM.R^2;
%plot((1:ind-2)*PARAM.checkpoint*PARAM.deltaT,risy(2*(n+m+j)+1,1:ind-2)/Uavg,'-')
plot((1:ind-2)*PARAM.checkpoint*PARAM.deltaT,risy(2*(n+m+j)+1,1:ind-2)/Umax,'-')
title('droplet velocity')
xlabel('time')
ylabel('v_{drop}/v_{max}')
grid on

%check force free
% subplot(2,2,3)
% Uavg = PARAM.Q/pi/PARAM.R^2;
% plot((1:ind-2)*PARAM.checkpoint*PARAM.deltaT,Fdrop(2:ind-1),'-')
% title('Force on drop')
% xlabel('time')
% ylabel('F')
% grid on

subplot(2,2,4)
V0 = 4/3*pi*PARAM.alpha^3*PARAM.R^3;
%plot volume variation
V_rel = abs(V-V0)/V0;
%figure
%semilogy((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
plot((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
xlabel('iteration')
ylabel('V_{rel}')
title('Volume variation')
grid on

if PARAM.continuation==1 || compareData==0
    disp('No comparison with Lac data')
elseif compareData==1
    figure
    ALPHA = [0.6 0.8 0.9 1 1.1 1.2 1.3];
    this = find(PARAM.alpha==ALPHA,1,'first');
    %compare with digitize
    name = '~/Documents/MATLAB/droplet_simulations/results/confined_droplet/validation/';
    filename = [name 'fig4a_lambda=' num2str(PARAM.visc) '.dat'];
    if PARAM.visc==0
        name = '~/Documents/MATLAB/droplet_simulations/results/confined_droplet/validationZeroLambda/';
        filename = [name 'alpha' num2str(PARAM.alpha) '_Lam0_z_and_r.txt'];
    end
    A = importdata(filename);
    xxx = A(:,1);
    yyy = A(:,2);
    if PARAM.visc==0
        xxx = xxx-max(xxx)+max(xDrop);
        plot([xxx; flip(xxx)],[yyy; -flip(yyy)])
        %hold on
        %plot(xxx,yyy)
    elseif PARAM.visc==1
        if this<5
            plot(xxx(1+(this-1)*101:101*this),yyy(1+(this-1)*101:101*this),'LineWidth',2)
        elseif this==5
            plot(xxx(402+(this-5)*201:402+201*(this-4)),yyy(402+(this-5)*201:402+201*(this-4)),'LineWidth',2)
        elseif this==6
            plot(xxx(602+(this-6)*101:602+101*(this-5)),yyy(602+(this-6)*101:602+101*(this-5)),'LineWidth',2)
        elseif this==7
            plot(xxx(702+(this-7)*101:702+201*(this-6)),yyy(702+(this-7)*201:702+201*(this-6)),'LineWidth',2)
        end
    elseif PARAM.visc==10
        if this<5
            plot(xxx(1+(this-1)*101:101*this),yyy(1+(this-1)*101:101*this),'LineWidth',2)
        elseif this==5
            plot(xxx(402+(this-5)*401:402+401*(this-4)),yyy(402+(this-5)*401:402+401*(this-4)),'LineWidth',2)
        elseif this==6
            plot(xxx(802+(this-6)*201:802+201*(this-5)),yyy(802+(this-6)*201:802+201*(this-5)),'LineWidth',2)
        elseif this==7
            plot(xxx(1002+(this-7)*101:1002+201*(this-6)),yyy(1002+(this-7)*201:1002+201*(this-6)),'LineWidth',2)
        end
    elseif PARAM.visc==0.1
        if this<6
            plot(xxx(1+(this-1)*101:101*this),yyy(1+(this-1)*101:101*this),'LineWidth',2)
        elseif this==6
            plot(xxx(502+(this-6)*201:502+201*(this-5)),yyy(502+(this-6)*201:502+201*(this-5)),'LineWidth',2)
        elseif this==7
            plot(xxx(702+(this-7)*101:702+201*(this-6)),yyy(702+(this-7)*201:702+201*(this-6)),'LineWidth',2)
        end
    end
    hold on
    plot(xDrop-max(xDrop),yDrop,'xr')
    plot(xDrop-max(xDrop),-yDrop,'xr')
    %plot(xWall,yWall,'k')
    %plot(xWall,-yWall,'k')
    axis([min(xDrop)-1 max(xDrop)+1 -2 2])
    axis equal
    %grid on
    title(['Validation Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ' \alpha=' num2str(PARAM.alpha)])
    hold off
    if PARAM.visc==0
        legend('Lailai','mine')
    else
        legend('Lac','mine')
    end
    grid on
end

