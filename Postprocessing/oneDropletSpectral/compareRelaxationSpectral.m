%post processing of test time stepping spectral

close all
clear variables

%set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%data
manyVisc = 1*ones(11,1);
%manyVisc = [0.01 0.03 0.1 0.3 1 3 10 31 100];
nnn = 50*ones(7,1);
manyTend = 100*ones(11,1);
%manyTend = [12 12 15 15 20 35 120 280 900];
ODE = 2;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1;         % 1 is extensional flow, 2 is rising droplet
volCorr = 0;
DT = 1e-2*ones(11,1);
%DT = [5e-4 5e-4 5e-4 5e-4 1e-3 1e-3 5e-3 1e-2 5e-2];
legendre = 1;
%manyD = [0.007 0.0162 0.0377 0.0874 0.2026 0.47];          % initial shape
manyD = [0.05 0.07 0.09 0.105:0.005:0.12];
%manyD = 0.091*ones(9,1);
CPUs = 16;
res = 1;        % results or server

if numel(nnn)==1
    cellLegend = cell(size(DT));
else
    cellLegend = cell(size(nnn));
end

%option
plotD = 1;
plotDelta = 1;
plotTimeScale = 1;
plotShape = 0;
semilogOn = 1;
zoom = 0;
plotModes = 0;
plotRes = 0;
plotAreaNorm = 0;
step = 1;
plotCurvSpace = 1;
plotPhaseSpacePolar = 0;    dimension = 2;  symmetric = 1;  modes = 4;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
end

%initialize
tScaleNum = zeros(numel(manyVisc),1);

for l = 1:numel(nnn)
    
    display([num2str(l) ' of ' num2str(numel(nnn))])
    
    %curent varibales
    n = nnn(l);
    visc = manyVisc(l);
    maxDT = DT(l);
    Tend = manyTend(l);
    Tbreak = 2*Tend;
    D = manyD(l);
    
    if nnn(1)~=nnn(2)
        %store element in cell for legend
        cellLegend(l) = {['n=' num2str(n)]};
    elseif DT(1)~=DT(2)
        %store element in cell for legend
        cellLegend(l) = {['dt=' num2str(n)]};
    elseif manyVisc(1)~=manyVisc(2)
        %store element in cell for legend
        cellLegend(l) = {['\lambda=' num2str(visc)]};
    elseif manyD(1)~=manyD(2)
        %store element in cell for legend
        cellLegend(l) = {['D=' num2str(D)]};
    end

        %filename
        name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=0_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(D) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '.mat'];

        %upload data
        upload = [dir name];
        load(upload)

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
        Dellipse = zeros(ite,1);
        delta = zeros(ite,1);
        Anorm = zeros(ite,1);
        V0 = 4/3*pi;
        k1_east = zeros(ite,1);
        k1_north = zeros(ite,1);
        k2_east = zeros(ite,1);
        k2_north = zeros(ite,1);

        %loops
        for i = 1:ite

           k = i*step;

           %display([num2str(i) ' of ' num2str(ite)])

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

           %from modes to grid
           [x,y] = fromModesToGrid(xMode,yMode,PARAM);
           
           %defromation parameter
           L = xInterp(0)-xInterp(1);
           B = 2*yInterp(0.5);
           Dellipse(i) = (L-B)/(L+B);
           delta(i) = L/B-1;

           %compute volume
           V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);

           R0 = nthroot(V(i)/4/pi*3,3);
           Anorm(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM) - 4*pi*R0^2;

           %compute force acting on droplet
           ForceFree(i) = ForceOnDropSpectralCurvilinear(x,y,PARAM.Ca,PARAM);
           
           if plotCurvSpace
               %compute curvature
               xp = PARAM.D1*x;    xpp = PARAM.D2*x;
               yp = PARAM.D1*y;    ypp = PARAM.D2*y;
               h = (xp.^2+yp.^2).^(0.5);
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
           end

           t = T;

           if plotShape==1

              xMid = (max(x)+min(x))/2;

              figure(1) 
              if BC==2
                    plot(y,-x,'kx-')
                    hold on
                    grid on
                    plot(-y,-x,'k')
                    title('Rising Droplet')
                    axis([-2 2 -2 2])
              elseif BC==1
                    plot(x,y,'kx-')
                    hold on
                    grid on
                    plot(x,-y,'k')
                    %axis([-2 2 -2 2])
                    title('Extensional')
              end
              xlabel('r')
              ylabel('x')
              axis equal
              %axis([-2 2 -2-xMid 2-xMid])
              hold off
              %title('Rising Droplet')

              if zoom==1

                 axis([-0.1 0.1 -1.2 -1.0]) 

              end

              drawnow

           end

           %stop when there are no more data
           if i~=ite
               if T(k+1)>Tbreak||(T(k+1)==0&&k>1)
                   ForceFree = ForceFree(1:i);
                   V = V(1:i);
                   t = T(1:i);
                   res = res(1:i);
                   Dellipse = Dellipse(1:i);
                   delta = delta(1:i);
                   Anorm = Anorm(1:i);
                   normElon = normElon(1:i);
                   k1_east = k1_east(1:i);
                   k1_north = k1_north(1:i);
                   k2_east = k2_east(1:i);
                   k2_north = k2_north(1:i);
                   disp('Break')
                   break;
               end
           end

        end
        
        %time scale
        Tr = (2*visc+3)*(19*visc+16)/(40*(1+visc));
        
        %get numerical time scale
        fitExp = fit(t(round(numel(t)/4):round(numel(t)/2)),Dellipse(round(numel(t)/4):round(numel(t)/2)),'exp1');
        tScaleNum(l) = -1/fitExp.b;
        
        if plotD==1

        %plot deformation parameter semilogy
        figure(2)
        if l>1
            hold on
        end
        Dnorm = Dellipse/Dellipse(1);
        plot(t/Tr,Dnorm)
        %xlabel('t/\tau_r')
        %ylabel('D/D_0')
        xyLabelTex('t/\tau_r','D/D_0')
        grid on
        axis([0 8 1e-4 1])
        title('Deformation parameter')
        
        if semilogOn==1
            
            figure(5)
            if l>1
                hold on
            end
            semilogy(t/Tr,Dellipse)
            %xlabel('t/\tau_r')
            %ylabel('D')
            xyLabelTex('t/\tau_r','D')
            grid on
            axis([0 8 1e-4 1])
            title('Deformation parameter')
            
        end
        
        end
        
        if plotAreaNorm==1

        %plot excess surface area semilogy
        figure(3)
        if l>1
            hold on
        end
        Anorm = Anorm/Anorm(1);
        plot(t,Dnorm)
        xlabel('t/\tau_r')
        ylabel('\Delta A/\DeltaA_0')
        grid on
        axis([0 8 0 1])
        title('Excess area')
        
        if semilogOn==1
            
            %plot excess surface area semilogy
            figure(6)
            if l>1
                hold on
            end
            Anorm = Anorm/Anorm(1);
            semilogy(t/Tr,Anorm)
            xlabel('t/\tau_r')
            ylabel('\Delta A/\DeltaA_0')
            grid on
            axis([0 8 0 1])
            title('Excess area')
            
        end
        
        end

        figure(4)
        if l>1
            hold on
        end
        errV = (V-V0)/V0;
        plot(t/Tr,errV)
        xlabel('t/\tau_r')
        ylabel('err_V')
        grid on
        title('Error on volume')
        
        %for ellispoidal relaxation
        [kNorthEllipse1,kEastEllipse1,kNorthEllipse2,kEastEllipse2] = ellipseCurvature(linspace(0,Dellipse(1),100));
        kNorthEllipse = kNorthEllipse1+kNorthEllipse2;
        kEastEllipse = kEastEllipse1+kEastEllipse2;
        k_east = k1_east+k2_east;
        k_north = k1_north+k2_north;
    
        %plot
        figure(7)
        plot(k1_north,k1_north+k1_east)
        hold on

        %plot
        figure(8)
        plot(k2_north,k2_north+k2_east)
        hold on
    
        %plot
        figure(9)
        plot(k_north,k_north+k_east)
        hold on
        
        if plotD==1

        %plot deformation parameter semilogy
        figure(10)
        if l>1
            hold on
        end
        deltaNorm = delta/delta(1);
        plot(t/Tr,deltaNorm)
        %xlabel('t/\tau_r')
        %ylabel('\delta/\delta_0')
        xyLabelTex('t/\tau_r','\delta/\delta_0')
        grid on
        axis([0 8 1e-4 1])
        title('Deformation parameter Campas')
        
        if semilogOn==1
            
            figure(11)
            if l>1
                hold on
            end
            semilogy(t/Tr,delta)
            %xlabel('t/\tau_r')
            %ylabel('\delta')
            xyLabelTex('t/\tau_r','\delta')
            grid on
            axis([0 8 1e-4 max(delta)])
            title('Deformation parameter Campas')
            
        end
        
        end

        %print simulation time
        disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])

    
end

if plotD==1
    figure(2)
    legend(cellLegend,'Location','Best')

    figure(5)
    legend(cellLegend,'Location','Best')
end

if plotDelta==1
    figure(10)
    legend(cellLegend,'Location','Best')

    figure(11)
    legend(cellLegend,'Location','Best')
end

if plotAreaNorm==1
    figure(3)
    legend(cellLegend,'Location','Best')

    figure(6)
    legend(cellLegend,'Location','Best')
end

figure(4)
legend(cellLegend,'Location','Best')

if plotCurvSpace==1
    
    figure(7)
    plot(kNorthEllipse1,kNorthEllipse1+kEastEllipse1,'k--')
    hold on
    plot(k1_north(1),k1_north(1)+k1_east(1),'.k','MarkerSize',50)
    plot(k1_north(end),k1_north(end)+k1_east(end),'k*','MarkerSize',20)
    %xlabel('k_1 north')
    %ylabel('k_1 north+k_1 east')
    xyLabelTex('K_1^{\rm{north}}','K_1^{\rm{north}}+K_1^{\rm{east}}')
    grid on
    cellLegend(l+1) = {'Ellipsoidal relaxation'};
    cellLegend(l+2) = {['t=' num2str(t(1))]};
    cellLegend(l+3) = {['t=' num2str(t(end))]};
    legend(cellLegend,'Location','Best')
    title('In plane curvature')
        
    figure(8)
    plot(kNorthEllipse2,kNorthEllipse2+kEastEllipse2,'k--')
    hold on
    plot(k2_north(1),k2_north(1)+k2_east(1),'.k','MarkerSize',50)
    plot(k2_north(end),k2_north(end)+k2_east(end),'k*','MarkerSize',20)
    %xlabel('k_2 north')
    %ylabel('k_2 north+k_2 east')
    xyLabelTex('K_2^{\rm{north}}','K_2^{\rm{north}}+K_2^{\rm{east}}')
    grid on
    cellLegend(l+1) = {'Ellipsoidal relaxation'};
    cellLegend(l+2) = {['t=' num2str(t(1))]};
    cellLegend(l+3) = {['t=' num2str(t(end))]};
    legend(cellLegend,'Location','Best')
    title('Azimuthal curvature')
    
    figure(9)
    plot(kNorthEllipse,kNorthEllipse+kEastEllipse,'k--')
    hold on
    plot(k_north(1),k_north(1)+k_east(1),'.k','MarkerSize',50)
    plot(k_north(end),k_north(end)+k_east(end),'k*','MarkerSize',20)
    %xlabel('k north')
    %ylabel('k north+k east')
    xyLabelTex('K^{\rm{north}}','K^{\rm{north}}+K^{\rm{east}}')
    grid on
    cellLegend(l+1) = {'Ellipsoidal relaxation'};
    cellLegend(l+2) = {['t=' num2str(t(1))]};
    cellLegend(l+3) = {['t=' num2str(t(end))]};
    legend(cellLegend,'Location','Best')
    title('Azimuthal curvature')
    title('Total curvature')
        
end

if plotTimeScale==1
    
    visc = logspace(-2,2,1000);
    Tr = (2*visc+3).*(19*visc+16)./(40*(1+visc));

    figure
    loglog(visc,Tr)
    xlabel('\lambda')
    ylabel('\tau/\tau_r')
    grid on
    hold on
    loglog(manyVisc,tScaleNum,'.','MarkerSize',40)
    legend('Analytical','Numerical','Location','Best')

end









