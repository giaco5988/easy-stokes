%post processing of bifurcation in extensional flow

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%parameters
n = 25;
lambda = 1;

%options
PlotLastShape = 0;
PlotRTheta = 0;
PlotCurv = 1;   polar = 0;
WhichShape = 50;

%load data
res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/';
%res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/server/';
%filename = [res 'continuation_extens_n=' num2str(n) '_visc=' num2str(lambda) '.mat'];
filename = [res 'ContSpectralExtensNewtonXY_n=' num2str(n) '_visc=' num2str(lambda) '.mat'];
load(filename)

%for which point I plot
%Dplot = [0.3 0.4 0.5 0.55 0.58 0.5834 0.595];
Dplot = 0.53:0.02:0.59;

ind0 = find(manyD==0,1);

manyD = manyD(1:ind0-1);
manyCa = manyCa(1:ind0-1);

%bifurcation diagram
figure
plot(manyCa,manyD,'k')
xlabel('Ca')
ylabel('D')
title('Bifurcation diagram')
grid on

for i = 1:numel(Dplot)
    
    display(i)
    
    [~,ind] = min(abs(Dplot(i)-manyD));
    x = manyA(:,ind);   y = manyB(:,ind);
    
    figure(1)
    hold on
    plot(manyCa(ind),manyD(ind),'.','MarkerSize',40)
    hold off
    
    figure(2)
    %subplot(2,3,1)
    hold on
    plot([x; flip(x)],[y; -flip(y)],'x-')
    grid on
    axis([-2.5 2.5 -1 1])
    axis equal
    hold off
    xlabel('x')
    ylabel('r')
    title('Shapes')
    
end

if PlotLastShape==1
    %plot last shape
    figure
    for i = 1:n+1

        %ind = find(manyA(1,:)==0,1,'first');
        x = manyA(:,WhichShape);   y = manyB(:,WhichShape);
        plot([0 x(i)],[0 y(i)])
        hold on
        axis equal
        grid on
        drawnow

    end
    plot(x,y)
    xlabel('x')
    ylabel('r')
    title(['Shape at iteration ' num2str(WhichShape)])
end

if PlotRTheta==1
    
    x = manyA(:,WhichShape);   y = manyB(:,WhichShape);
    
    %plot last r(theta)
    figure
    plot(theta,sqrt(x.^2+y.^2),'o-')
    grid on
    xlabel('\theta')
    ylabel('R(\theta)')
    title(['Solution at iteration ' num2str(WhichShape)])
end

if PlotCurv==1
    
    %coordinates
    x = manyA(:,WhichShape);   y = manyB(:,WhichShape);
    
    if polar==1
    
        %compute geometry
        r = sqrt(x.^2+y.^2);
        %compute rho in symmetry axis
        fVolume = @(rho) ModifyVolumeSpectral(theta,r,rho,V0,Replace,thetaReplace,PARAM);
        options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
        rMiddle = fsolve(fVolume,1,options);

        r = [r(1:Replace-1); rMiddle; r(Replace:end)];
        theta = [theta(1:Replace-1); thetaReplace; theta(Replace:end)];

        D1 = PARAM.D1;
        D2 = PARAM.D2;
        rp = D1*r;    rpp = D2*r;

        %compute curvature in meridional plane
        h = sqrt(r.^2+rp.^2);
        K1 = (r.^2+2*rp.^2-r.*rpp)./h.^3;
    
    elseif polar==0
        
        %compute rho in symmetry axis
        D1 = PARAM.D1;
        D2 = PARAM.D2;
        
        %compute geomtrical derivaties
        xp = D1*x;    yp = D1*y;
        xpp = D2*x;    ypp = D2*y;

        %compute normal vector
        h = (xp.^2+yp.^2).^(0.5);
        nx = yp./h;
        ny = -xp./h;
        
        %compute curvature in meridional plane
        K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
        
    end
    
    %plot last r(theta)
    figure
    plot(x,K1,'o-')
    grid on
    xlabel('\theta')
    ylabel('K_1')
    title(['Curvature at iteration ' num2str(WhichShape)])
end

