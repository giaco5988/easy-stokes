%post processing of test time stepping spectral

close all
clear variables

%data
Ca = 5;
visc = 1;
n = 50;
Tend = 20;
ODE = 2;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 2;         % 1 is extensional flow, 2 is rising droplet
volCorr = 0;
maxDT = 1e-3;
legendre = 1;
D = [0.2 0.1 0.3];          % initial shape
CPUs = 16;
res = 1;        %results, server or relaxationCampas

%option
plotShape = 1;  threeD = 0;
plotSnapshot = 0;   timeSnapshot = [0 10 20 26];
savePlotShape = 0;  nameSave = 'extensNoBreak';     dest = '/Users/Giacomo/Documents/research_notes/presentationPHDdefense/movie/frame/';
getFewShapes = 0;   Tget = linspace(0,18,5);
plotOnlyNorm = 1;
plotOnlyArea = 0;
zoom = 0;
plotModes = 0;
plotRes = 0;    plotPower = 1;
plotAreaNorm = 1;
Tbreak = 2*Tend;
step = 1;
plotPhaseSpacePolar = 0;    dimension = 2;  symmetric = 1;  modes = 4;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==2
    dir = '~/Documents/MATLAB/droplet_simulations/results/relaxationCampas/';
end

%filename
%name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '.mat'];
name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(D) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '.mat'];
%name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(D) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '_CPUs=' num2str(CPUs) '.mat'];

%upload data
upload = [dir name];
load(upload)

%PARAM.legendre = 0;

%number of iteration
ite = round(numel(T)/step);

%initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%initialize
ForceFree = zeros(ite,1);
normElon = zeros(ite,1);
res = zeros(ite,1);
V = zeros(ite,1);
Area = zeros(ite,1);
Evisc = zeros(ite,1);
Dellipse = zeros(ite,1);
xcm = zeros(ite,1);
Anorm = zeros(ite,1);
k1_east = zeros(ite,1);
k1_north = zeros(ite,1);
k2_east = zeros(ite,1);
k2_north = zeros(ite,1);
f1 = zeros(ite,1);
f2 = zeros(ite,1);
f3 = zeros(ite,1);
V0 = 4/3*pi;
xGetShape = cell(size(Tget));
yGetShape = cell(size(Tget));
countGet = 1;
indDots = zeros(numel(timeSnapshot),1);

%loops
plotCount = 1;
shift = 0;  islast = 0;
for i = 1:ite
    
   k = i*step;
    
   display([num2str(i) ' of ' num2str(ite)])
    
   %get currents modes
   xyMode = Y(k,:)';
   xMode = Y(k,1:2:end-1)';
   yMode = Y(k,2:2:end)';
   
   %interpolation
    if PARAM.legendre==2 || PARAM.legendre==1

          xInterp = @(t) interpLegendreZeroOne(t,xMode);
          yInterp = @(t) interpLegendreZeroOne(t,yMode);

    elseif PARAM.legendre==0

          xInterp = chebcoeffs2chebvals(xMode);
          yInterp = chebcoeffs2chebvals(yMode);
          
          xInterp = chebfun(xInterp,[0 1]);
          yInterp = chebfun(yInterp,[0 1]);
          
    end
   
   [x,y] = fromModesToGrid(xMode,yMode,PARAM);
   
   if plotPhaseSpacePolar==1
            
            %radius modes
            f = LegendreSerie(x,y,modes,symmetric,PARAM);
            f1(i) = f(2);
            f2(i) = f(3);
            f3(i) = f(4);
        
   end
   
   %compute elongation norm
   RforNorm = sqrt(x.^2+y.^2);
   xcm(i) = CenterMassCurvAxisSpectral(x,y,PARAM);
   if Ca>=0
        normElon(i) = abs(RforNorm(1)-xcm(i));
   else
       [~,ind] = max(abs(y));
        normElon(i) = y(ind);
   end
   
   %compute volume abd surface area
   V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);
   Area(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM);
   
   if plotAreaNorm==1
       
       R0 = nthroot(V(i)/4/pi*3,3);
       Anorm(i) = (Area(i) - 4*pi*R0^2)./(Area(1) - 4*pi*R0^2);
       
   end
   
   if plotSnapshot==1
   
      if islast==0
          if timeSnapshot(plotCount)==T(i)

              indDots(plotCount) = i;
              totalShapes = numel(timeSnapshot);

              figure(10)
              plot(x,y+shift,'k')
              hold on
              plot(x,-y+shift,'k')
              %grid on
              axis equal
              xlabel('z')
              ylabel('r')
              axis([-5 5 -2 2+3*(totalShapes-1)])

              plotCount = plotCount+1;
              shift = shift+3;
              axis off

              if totalShapes==plotCount-1
                  islast = 1;
              end

          end
      end
   
   end
   
   %compute force acting on droplet
   ForceFree(i) = ForceOnDropSpectralCurvilinear(x,y,PARAM.Ca,PARAM);
   
   %defromation parameter
   %L = max(x)-min(x);
   L = xInterp(0)-xInterp(1);
   %B = 2*y((PARAM.n+1)/2);
   B = 2*yInterp(0.5);
   Dellipse(i) = (L-B)/(L+B);
   
   %compute curvature
   xp = PARAM.D1*x;    xpp = PARAM.D2*x;
   yp = PARAM.D1*y;    ypp = PARAM.D2*y;
   h = (xp.^2+yp.^2).^(0.5);
   nx = yp./h;
   ny = -xp./h;
   K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
   if PARAM.legendre==1
        K2 = ny./y;
   elseif PARAM.legendre==0||PARAM.legendre==2
        K2 = ny./y;
        K2([1 end]) = K1([1 end]);
   end
   
   %curvature on tip and top
   k1_east(i) = K1(1);
   k1_north(i) = K1(round(PARAM.n/2));
   k2_east(i) = K2(1);
   k2_north(i) = K2(round(PARAM.n/2));
   
   %compute viscous dissipation
   fx = (K1+K2).*nx;    fy = (K1+K2).*ny;
   
   
   %get time
   t = T;
   
   if getFewShapes==1
       
       if sum(t(i)==Tget)
           
            xGetShape(countGet) = {x};
            yGetShape(countGet) = {y};
            
            countGet = countGet+1;
           
       end
       
   end
   
   if plotShape==1
       
      xMid = (max(x)+min(x))/2;
      
      if savePlotShape==1 && i==1
        width = 1201;
        height = 500;
        figure('pos',[50 100 width height])
      end
      
      figure(1)
      if threeD==1
        
          solidOfRevolution(y',x',[0 2*pi],[0 0.8 1],0.8)
          axis off
          axis equal
          axis([-5 5 -1 1 -1 1])
        
      else
          if BC==2
                plot(y,-x,'k-')
                hold on
                grid on
                plot(-y,-x,'k')
                title('Rising Droplet')
                axis([-2 2 -2 2])
          elseif BC==1 || BC==3
                plot(x,y,'k-')
                hold on
                grid on
                plot(x,-y,'k')
                axis([-4 4 -2 2])
                title('Extensional')
          end
          xlabel('r')
          ylabel('x')
          axis equal
          %axis([-2 2 -2-xMid 2-xMid])
          hold off
          %title('Rising Droplet')
      end
      
      if savePlotShape==1
          
        disp('Save')
        grid off
        axis off
        title('')
        print('-dpng','-loose','-r400',[dest nameSave sprintf('%03d',i-1) '.png'])
          
      end
      
      if zoom==1
          
         axis([-0.1 0.1 -1.2 -1.0]) 
          
      end
      
      drawnow
      
   end
   
   if plotModes==1
        %plot modes
        figure(2)
        loglog(abs(xMode(2:2:end)))
        hold on
        loglog(abs(yMode(1:2:end)))
        grid on
        hold off
        xlabel('n')
        ylabel('c')
        legend('C_x','C_y')
        drawnow
   end
   
   if plotRes==1 || plotPower==1
       
           
           if PARAM.legendre==1||PARAM.legendre==2
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [uvModes,res(i)] = dropLegendreCurvilinearModes(t,xyMode,PARAM);
               cd(here)
           elseif PARAM.legendre==0
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [uvModes,res(i)] = dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
               cd(here)
           end
           
           [u,v] = fromModesToGrid(uvModes(1:2:end-1),uvModes(2:2:end),PARAM);
           fE = fx.*u + fy.*v;
           Evisc(i) = axisIntSpectralCurvilinear(x,y,fE,PARAM);
       
   end
   
   %stop when there are no more data
   if i~=ite
   if T(k+1)>Tbreak||(T(k+1)==0&&k>1)
       ForceFree = ForceFree(1:i);
       V = V(1:i);
       t = T(1:i);
       res = res(1:i);
       normElon = normElon(1:i);
       disp('Break')
       break;
   end
   end
    
end

figure
%subplot(3,2,1)
if BC==2
    plot(t,normElon,'k-')
    grid on
    title('Rising Droplet')
    xlabel('t')
    ylabel('||r(0)-1||_{\infty}')
elseif BC==1
    plot(t,normElon,'k-')
    grid on
    title('Extensional flow')
    xlabel('t')
    ylabel('||r(0)-1||_{\infty}')
end

%subplot(3,2,2)
figure
loglog(abs(xMode(2:2:end)))
hold on
loglog(abs(yMode(1:2:end)))
grid on
hold off
xlabel('n')
ylabel('c')
legend('C_x','C_y')

errV = (V-V0)/V0;
%subplot(3,2,3)
figure
plot(t,errV)
xlabel('t')
ylabel('err_V')
grid on
title('Error On Volume')

% subplot(3,2,4)
% plot(t,ForceFree)
% xlabel('t')
% ylabel('F')
% title('Force On Drop')
% grid on

%subplot(3,2,4)
figure
plot(t(1:i),Anorm(1:i))
xlabel('t')
ylabel('\Delta A')
title('Excess area')
grid on

%subplot(3,2,5)
figure(5)
semilogy(t,res)
hold on
xlabel('t')
ylabel('res')
title('residuals')
grid on

%compute curvature
xp = PARAM.D1*x;    xpp = PARAM.D2*x;
yp = PARAM.D1*y;    ypp = PARAM.D2*y;
h = (xp.^2+yp.^2).^(0.5);
nx = yp./h;
ny = -xp./h;
K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
if PARAM.legendre==1
    K2 = ny./y;
elseif PARAM.legendre==0||PARAM.legendre==2
    K2 = ny./y;
    K2([1 end]) = K1([1 end]);
end

%subplot(3,2,6)
figure
plot(x,K1)
hold on
plot(x,K2)
plot(x,K1+K2,'k')
xlabel('x')
ylabel('k')
legend('k_1','k_2','k','Location','Best')
title('curvature')
grid on

%plot for droplet relaxation
if Ca==0
    
    %time scale
    Tr = (2*visc+3)*(19*visc+16)/(40*(1+visc));
    
    %plot deformation parameter
    figure
    %subplot(2,1,1)
    Dnorm = Dellipse/Dellipse(1);
    plot(t/Tr,Dnorm(1:numel(t)),'k')
    xlabel('t/\tau_r')
    ylabel('D/D_0')
    grid on
    title(['Drop relaxation \lambda=' num2str(visc)])
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          plot(t(ind)/Tr,Dnorm(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot deformation parameter semilogy
    figure
    %subplot(2,1,2)
    Dnorm = Dellipse/Dellipse(1);
    semilogy(t/Tr,Dnorm(1:numel(t)),'k')
    xlabel('t/\tau_r')
    ylabel('D/D_0')
    grid on
    title(['Drop relaxation \lambda=' num2str(visc)])
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          semilogy(t(ind)/Tr,Dnorm(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot excess surface area
    figure
    %subplot(2,1,1)
    Anorm = Anorm/Anorm(1);
    plot(t/Tr,Anorm(1:numel(t)),'k')
    xlabel('t/\tau_r')
    ylabel('\Delta A/\DeltaA_0')
    grid on
    title(['Drop relaxation \lambda=' num2str(visc)])
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          plot(t(ind)/Tr,Anorm(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot excess surface area semilogy
    figure
    %subplot(2,1,2)
    Anorm = Anorm/Anorm(1);
    semilogy(t/Tr,Anorm(1:numel(t)),'k')
    xlabel('t/\tau_r')
    ylabel('\Delta A/\DeltaA_0')
    grid on
    title(['Drop relaxation \lambda=' num2str(visc)])
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          semilogy(t(ind)/Tr,Anorm(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %for ellispoidal relaxation
    [kNorthEllipse1,kEastEllipse1,kNorthEllipse2,kEastEllipse2] = ellipseCurvature(linspace(0,Dellipse(1),100));
    kNorthEllipse = kNorthEllipse1+kNorthEllipse2;
    kEastEllipse = kEastEllipse1+kEastEllipse2;
    k_east = k1_east+k2_east;
    k_north = k1_north+k2_north;
    
    %plot
    figure
    plot(k1_north,k1_north+k1_east,'k')
    hold on
    plot(kNorthEllipse1,kNorthEllipse1+kEastEllipse1,'k--')
%     plot(k1_north(1),k1_north(1)+k1_east(1),'.','MarkerSize',40)
%     plot(k1_north(end),k1_north(end)+k1_east(end),'.','MarkerSize',40)
    xlabel('k_1 north')
    ylabel('k_1 north+k_1 east')
    grid on
    %legend('Data from simulation','Ellipsoidal relaxation','t=0',['t=' num2str(t(end))],'Location','Best')
    legend('Data from simulation','Ellipsoidal relaxation','Location','Best')
    title('In plane curvature')
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          plot(k1_north(ind),k1_north(ind)+k1_east(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot
    figure
    plot(k2_north,k2_north+k2_east,'k')
    hold on
    plot(kNorthEllipse2,kNorthEllipse2+kEastEllipse2,'k--')
%     plot(k2_north(1),k2_north(1)+k2_east(1),'.','MarkerSize',40)
%     plot(k2_north(end),k2_north(end)+k2_east(end),'.','MarkerSize',40)
    xlabel('k_2 north')
    ylabel('k_2 north+k_2 east')
    grid on
    %legend('Data from simulation','Ellipsoidal relaxation','t=0',['t=' num2str(t(end))],'Location','Best')
    legend('Data from simulation','Ellipsoidal relaxation','Location','Best')
    title('Azimuthal curvature')
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          plot(k2_north(ind),k2_north(ind)+k2_east(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot
    figure
    plot(k_north,k_north+k_east,'k')
    hold on
    plot(kNorthEllipse,kNorthEllipse+kEastEllipse,'k--')
%     plot(k_north(1),k_north(1)+k_east(1),'.','MarkerSize',40)
%     plot(k_north(end),k_north(end)+k_east(end),'.','MarkerSize',40)
    xlabel('k north')
    ylabel('k north+k east')
    grid on
    %legend('Data from simulation','Ellipsoidal relaxation','t=0',['t=' num2str(t(end))],'Location','Best')
    legend('Data from simulation','Ellipsoidal relaxation','Location','Best')
    title('Azimuthal curvature')
    title('Total curvature')
    
    if getFewShapes==1
        
       hold on
       for i = 1:numel(Tget)
           
          ind = find(Tget(i)==t,1,'first');
          plot(k_north(ind),k_north(ind)+k_east(ind),'.','MarkerSize',40)
           
       end
        
    end
    
    %plot few shapes
    figure
    shift = 0;
    for i = 1:numel(Tget)
           
        ind = find(Tget(i)==t,1,'first');
        hold on
        x = xGetShape{i};   y = yGetShape{i};
        plot([x; flip(x)],[y; -flip(y)]-shift)
        axis equal
        grid on
        xlabel('x')
        ylabel('r')
        title('Droplet relaxation')
        axis([-3.5 3.5 -1.5 1.5])
        %axis off
        
        %shift = shift+3;
           
    end
    
    if plotPower==1
    
        %compute energy balance
        D1 = finiteDifference1D(numel(t),[2 0],1);
        dAdt = D1*Area(1:numel(t))./(D1*t);

        figure
        plot(t,-dAdt)
        hold on
        grid on
        plot(t,-Evisc(1:numel(t)),'--')
        xlabel('t')
        ylabel('E_s,-E_v')
        title('Power balance')
        legend('Surface area','Viscous dissip','Location','Best')
    
    end
    
end

%plot phase space
if plotPhaseSpacePolar==1
        
        figure(11)
        
        if dimension==2
            
            hold on
            plot(f1(1:i),f2(1:i))
            grid on
            plot(f1(indDots),f2(indDots),'k.','MarkerSize',30)
            
            
        elseif dimension==3
            
            hold on
            plot3(f1(1:i),f2(1:i),f3(1:i))
            zlabel('f_3')
            grid on
            
        end
        
        xlabel('f_1')
        ylabel('f_2')
        title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
end

if plotOnlyNorm==1

    figure
    %semilogy(t,normElon,'k')
    plot(t,normElon,'k')
    xlabel('t')
    ylabel('L/a')
    title('Stone elongation')
    grid on
    %axis([0 200 1 10])
    
    if plotSnapshot==1
        
        hold on
        plot(t(indDots),normElon(indDots),'k.','MarkerSize',30)
        
        
    end
    
end

if plotOnlyArea==1

    figure
    semilogy(t(1:i),Anorm(1:i),'k')
    xlabel('t')
    ylabel('\Delta S')
    title('Excess area')
    grid on

end

%print simulation time
disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











