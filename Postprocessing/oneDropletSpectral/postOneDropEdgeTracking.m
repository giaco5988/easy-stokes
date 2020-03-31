%post processing of test time stepping spectral

close all
clear variables

bemDir = '~/Documents/MATLAB/droplet_simulations/dropSpectral';

%data
Ca = 0.25;
visc = 1;
n = 50;
ODE = 2;            % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 2;             % 1 is extensional flow, 2 is rising droplet
volCorr = 1;
maxDT = 2e-4;
legendre = 1;
CPUs = 16;
res = 0;            % 0 server, 1 results, 2 is edge tracking dir, 3 is paper dir
edgeLoop = 1000;
edgeDelta = 1e-4;

%option
rangeLoop = 501:501;
%rangeLoop = [11 12 22 23];   %this for paper 2 
%rangeLoop = [11 12];
%rangeLoop = [1 2 15 14 25 26];
%rangeLoop = [1 2 15 14];
%rangeLoop = [1 2 12 13]; %this for paper
%rangeLoop = [1 2 8 6];
colIN = '--'; colOUT = '-';
plotShape = 1;   mesh = 1;  saveShape = 0;
zoomShape = 0;   meshZoom = 0;
plotNorm = 1;
plotL2Norm = 0;
plotModes = 0;
plotRes = 0;
plotDefParameter = 0;
plotCurv = 0;
plotVol = 0;
plotExcessArea = 0;
plotForce = 0;
plotPhaseSpacePolar = 0;    dimension = 2;  symmetric = 0;  modes = 7;
convergeShape  = 0.02;

%snapshot
plotSnapshots = 0;
stopPlot = 1;   stopPlotEnd = [750 500];
% cellShapes{1} = [0:20:40 75];       %paper
% cellShapes{2} = [0:20:40 45];       %paper
% cellShapes{3} = 0:30:150;           %paper
% cellShapes{4} = [1.7:31:95 120];    %paper
cellShapes{1} = [0 30 60 90 115];       %paper 2
cellShapes{2} = [0 30 60 90 120 150];       %paper 2
cellShapes{3} = [68.8 100 130 180 205];           %paper 2
cellShapes{4} = [75.8 100 130 180 210 240];    %paper 2

%saving options
plotMovie = 0;
saveData = [0 0 0 0];
plotLoop = 2;
stepSave = 10;
savePlotEnd = [115 180 205 240];
solidLine = [1 1 1 1];
nameShape = 'edgeShapeTwo';
nameNorm = 'elongNormTwo';
nameRes = 'residuaslTwo';
namePhase = 'phaseSpaceTwo';

%plot frequency
Tbreak = Inf;
step = 1;

%destination for saving plot
saveDest = '/Users/Giacomo/Documents/research_notes/EuroMecBari/presentation/movie/frame/';

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==2
    dir = '~/Documents/MATLAB/droplet_simulations/results/edgeTracking/';
elseif res==3
    dir = '~/Documents/MATLAB/droplet_simulations/results/edgeTracking/forPaper/data/';
end

%filename
name = ['edgeTrackingDrop_edgeLoop=' num2str(edgeLoop) '_deltaEdge=' num2str(edgeDelta) '_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(volCorr) '.mat'];

%upload data
upload = [dir name];
load(upload)

%initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%iintialize
if plotMovie==1
    
   Tstore = cell(numel(rangeLoop),1);
   Ystore = cell(numel(rangeLoop),1);
   ElongStore = cell(numel(rangeLoop),1);
   f1Store = cell(numel(rangeLoop),1);
   f2Store = cell(numel(rangeLoop),1);
   resStore = cell(numel(rangeLoop),1);
    
end

%get colot rgb
color = get(gca,'ColorOrder');
close gcf
countShapeLoop = 1;
for l = rangeLoop
    
    
    T = Tedge{l};
    Y = Yedge{l};
    if plotSnapshots==1
        Tshape = cellShapes{countShapeLoop};
        indDots = zeros(numel(Tshape),1);
    end
    countShapeLoop = countShapeLoop+1;
    
    display([num2str(l) ' of ' num2str(numel(Tedge))])
    
    if isempty(T)
        warning('No data')
        break;
    end

    %number of iteration
    ite = round(numel(T)/step);

    %initialize
    ForceFree = zeros(ite,1);
    normElon = zeros(ite,1);
    normL2 = zeros(ite,1);
    res = zeros(ite,1);
    V = zeros(ite,1);
    DDD = zeros(ite,1);
    A = zeros(ite,1);
    V0 = 4/3*pi;
    f1 = zeros(ite,1);
    f2 = zeros(ite,1);
    f3 = zeros(ite,1);

    %loops
    countShift = 1;
    countShape = 1;
    for i = 1:ite

       k = i*step;

       %get currents modes
       xyMode = Y(k,:)';
       xMode = Y(k,1:2:end-1)';
       yMode = Y(k,2:2:end)';

       if PARAM.legendre==1||PARAM.legendre==2

           x = LegendreBuildXY(xMode,PARAM.PPP);
           y = LegendreBuildXY(yMode,PARAM.PPP);

       elseif PARAM.legendre==0

           x = chebcoeffs2chebvals(xMode);
           y = chebcoeffs2chebvals(yMode);

       end
       
       if plotPhaseSpacePolar==1
            
            %radius modes
            f = LegendreSerie(x,y,modes,symmetric,PARAM);
            f1(i) = f(3);
            f2(i) = f(5);
            f3(i) = f(7);
        
       end

       %compute volume
       V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);
       xcm = CenterMassCurvAxisSpectral(x,y,PARAM);
       R0 = nthroot(3/4*V(i)/pi,3);
       
       %deformation parameter
       L = max(x)-min(x);
       B = 2*y(round(numel(y)/2));
       DDD(k) = (L-B)/(L+B);
       
       %compute elongation norm
       RforNorm = sqrt(x.^2+y.^2);
       normElon(i) = x(1)-xcm;
       %normElon(i) = DDD(k);
       normL2(i) = sqrt(PARAM.WG'*(RforNorm-1).^2);
       
       %compute area
       A(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM);

       %compute force acting on droplet
       ForceFree(i) = ForceOnDropSpectralCurvilinear(x,y,PARAM.Ca,PARAM);

       t = T;
       
       if plotSnapshots==1 && sum(abs(T(i)-Tshape)<1e-4)
              
              shiftX = 15;
              shiftY = 2;
              
              figure(12)
              plot([x; flip(x)]+shiftX*(countShapeLoop-1),[y; -flip(y)]+shiftY*(countShift-1),'Color',color(countShapeLoop-1,:))
              %plot(x,y,'Color',color(l,:))
              axis off
              axis equal
              hold on
              title('Shapes')
              
              %find index for dots
              indDots(countShape) = i;
              
              %counters
              countShift = countShift+1;
              countShape = countShape+1;
               
       end

       if plotShape==1

          xMid = (max(x)+min(x))/2;

          figure(1) 
          if BC==2
                plot(y,-x,'k-')
                hold on
                grid on
                plot(-y,-x,'k')
                title('Rising Droplet')
                axis([-2 2 -2-xcm 2-xcm])
                if mesh==1
                  plot(y,-x,'xr')
                end
          elseif BC==1
                plot(x,y,'k-')
                hold on
                grid on
                plot(x,-y,'k')
                axis([-3 3 -3 3])
                title(['Edge shape Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
                if mesh==1
                  plot(x,-y,'xr')
                end
          end
          xlabel('r')
          ylabel('x')
          axis equal
          if zoomShape==0
            hold off
            drawnow
          end
          %title('Rising Droplet')
          
          if saveShape==1 && countShapeLoop-1==saveLoop
          
            hold off
            plot(x,y,'Color',color(countShapeLoop-1,:))
            hold on
            plot(x,-y,'Color',color(countShapeLoop-1,:))
            disp('Save')
            axis equal
            axis([-10 10 -2 2])
            grid off
            axis off
            title('')
            print('-dpng','-loose','-r100',[saveDest nameShape sprintf('%03d',i-2) '.png'])

          end
          
       end
       
       if zoomShape==1
           
           if plotShape==1
             figure(1)
             plot([-0.1 0.1],[-1.15 -1.15],'r')
             plot([0.1 0.1],[-1.15 -0.95],'r')
             plot([-0.1 0.1],[-0.95 -0.95],'r')
             plot([-0.1 -0.1],[-1.15 -0.95],'r')
             hold off
             drawnow
           end
             
             figure(11)
             if BC==2
                plot(y,-x,'k-')
                hold on
                grid on
                plot(-y,-x,'k')
                title('Rising Droplet')
                axis([-2 2 -2 2])
                if meshZoom==1
                  plot(y,-x,'xr')
                end
                xlabel('r')
                ylabel('x')
                axis equal
                hold off
                axis([-0.1 0.1 -1.15 -0.95])
             
             end

           drawnow

       end

       if plotRes==1

               here = pwd;
               cd(bemDir)
               if PARAM.legendre==1
                   [~,res(i)] = dropLegendreCurvilinearModes(t,xyMode,PARAM);
               elseif PARAM.legendre==0
                   [~,res(i)] = dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
               end
               cd(here)

       end

       %stop when there are no more data
       if i~=ite
           if T(k+1)>Tbreak||(T(k+1)==0&&k>1)
               ForceFree = ForceFree(1:i);
               V = V(1:i);
               t = T(1:i);
               res = res(1:i);
               normElon = normElon(1:i);
               normL2 = normL2(1:i);
               disp('Break')
               break;
           end
       end

    end
    
    if stopPlot==1 && plotSnapshots==1
        %finalInd = find(abs(stopPlotTime(countShapeLoop)-t),1,'first');
        finalInd = indDots(end); 
    else
        finalInd = numel(t);
    end
    
    %figure out if it is in or out the bassin of attraction
    L = max(x)-min(x);
    B = 2*y(round(numel(y)/2));
    D = (L-B)/(L+B);
    
    if PARAM.BC==1
    
        here = pwd;
        cd(bemDir)
        if PARAM.visc==1
                    load('./steadyState/CaExt')
                    load('./steadyState/DExt')
        elseif PARAM.visc==0
                    load('./steadyState/CaExt0')
                    load('./steadyState/DExt0')
        elseif PARAM.visc==0.02
                    load('./steadyState/CaExt002')
                    load('./steadyState/DExt002')
        elseif PARAM.visc==0.05
                    load('./steadyState/CaExt005')
                    load('./steadyState/DExt005')
        elseif PARAM.visc==0.1
                    load('./steadyState/CaExt01')
                    load('./steadyState/DExt01')
        elseif PARAM.visc==0.5
                    load('./steadyState/CaExt05')
                    load('./steadyState/DExt05')
        elseif PARAM.visc==5
                    load('./steadyState/CaExt5')
                    load('./steadyState/DExt5')
        elseif PARAM.visc==10
                    load('./steadyState/CaExt10')
                    load('./steadyState/DExt10')
        elseif PARAM.visc==0.01
                    load('./steadyState/CaExt001')
                    load('./steadyState/DExt001')
        end
        cd(here)
        [~,ind] = min(abs(PARAM.Ca-manyCa));
        Dca = manyD(ind);
        OutIn  = abs(D-Dca)/Dca > convergeShape;
    
    elseif PARAM.BC==2

        OutIn  = abs(D) > convergeShape;

    end
    
    if OutIn==1
        col = colOUT;
    elseif OutIn==0
        col = colIN;
    end

    if plotModes==1
        figure(4)
        loglog(abs(xMode(2:2:end)))
        hold on
        loglog(abs(yMode(1:2:end)))
        grid on
        hold off
        xlabel('n')
        ylabel('c')
        legend('C_x','C_y')
    end
    
    if plotVol==1
        errV = (V-V0)/V0;
        figure(5)
        plot(t,errV)
        hold on
        xlabel('t')
        ylabel('err_V')
        grid on
        title('Error On Volume')
        drawnow
    end
    
    if plotExcessArea==1
        excess = A-4*pi*R0^2;
        %excess = A-4*pi;
        %excess = A;
        figure(10)
        plot(t,excess,col)
        hold on
        xlabel('t')
        ylabel('A-4\pi')
        grid on
        title('Ecxess area')
        drawnow
    end
    
    if plotForce==1
        figure(6)
        plot(t,ForceFree)
        hold on
        xlabel('t')
        ylabel('F')
        title('Force On Drop')
        grid on
    end
    
    if plotRes==1
        
        figure(20)
        semilogy(t,res,col)
        hold on
        xlabel('t')
        ylabel('res')
        title(['residuals Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        grid on
        drawnow
        
        if plotSnapshots
                
               hold on
               semilogy(t(indDots),res(indDots),'.k','MarkerSize',22) 
                
        end
        
    end

    if plotCurv==1
        
        %compute curvature
        xp = PARAM.D1*x;    xpp = PARAM.D2*x;
        yp = PARAM.D1*y;    ypp = PARAM.D2*y;
        h = (xp.^2+yp.^2).^(0.5);
        nx = yp./h;
        ny = -xp./h;
        K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
        if PARAM.legendre==1
            K2 = ny./y;
        elseif PARAM.legendre==0 || PARAM.legendre==2
            K2 = ny./y;
            K2([1 end]) = K1([1 end]);
        end
    
        figure(8)
        plot(x,K1+K2)
        %hold on
        %plot(x,K2)
        %xlabel('x')
        %ylabel('k')
        xyLabelTex('z','K')
        %legend('k_1','k_2')
        title('curvature')
        grid on
        hold off
    end
    
    if plotL2Norm==1
        figure(9)
        if BC==2
            plot(t,normL2,col)
            hold on
            grid on
            title('Rising Droplet')
            xlabel('t')
            ylabel('||r(0)-1||_2')
        elseif BC==1
            plot(t,normL2,col)
            grid on
            hold on
            title('Extensional flow')
            xlabel('t')
            ylabel('||r-1||_2')
        end
        drawnow
    end
    
    if plotNorm==1
        figure(3)
            
        if BC==2
            semilogy(t(1:finalInd),normElon(1:finalInd),col)
            hold on
            grid on
            title('Rising Droplet')
            xlabel('t')
            ylabel('L/a')
        elseif BC==1
            semilogy(t(1:finalInd),normElon(1:finalInd),col)
            %semilogy(t,normElon,col)
            grid on
            hold on
            title(['Extensional flow Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
            xlabel('t')
            ylabel('L/a')
        end
        drawnow
        
        if plotSnapshots==1
                
               semilogy(t(indDots),normElon(indDots),'.k','MarkerSize',22) 
                
        end
    end
    
    %plot phase space
    if plotPhaseSpacePolar==1
        
        figure(11)
        
        if dimension==2
            
            if l>rangeLoop(1)
                hold on
            end
            plot(f1(1:finalInd),f2(1:finalInd),col)
            grid on
            
            if plotSnapshots==1
                
               plot(f1(indDots),f2(indDots),'.k','MarkerSize',22) 
                
            end
            
        elseif dimension==3
            
            hold on
            plot3(f1(1:finalInd),f2(1:finalInd),f3(1:finalInd),col)
            zlabel('f_3')
            grid on
            
        end
        
        xlabel('f_2')
        ylabel('f_4')
        title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
    end
    
    %plot phase space
    if plotDefParameter==1
        
        if PARAM.BC==1
            here = pwd;
            cd(bemDir)
            if PARAM.visc==1
                    load('./steadyState/CaExt')
                    load('./steadyState/DExt')
            elseif PARAM.visc==0
                    load('./steadyState/CaExt0')
                    load('./steadyState/DExt0')
            elseif PARAM.visc==0.02
                    load('./steadyState/CaExt002')
                    load('./steadyState/DExt002')
            elseif PARAM.visc==0.05
                    load('./steadyState/CaExt005')
                    load('./steadyState/DExt005')
            elseif PARAM.visc==0.1
                    load('./steadyState/CaExt01')
                    load('./steadyState/DExt01')
            elseif PARAM.visc==0.5
                    load('./steadyState/CaExt05')
                    load('./steadyState/DExt05')
            elseif PARAM.visc==5
                    load('./steadyState/CaExt5')
                    load('./steadyState/DExt5')
            elseif PARAM.visc==10
                    load('./steadyState/CaExt10')
                    load('./steadyState/DExt10')
            elseif PARAM.visc==0.01
                    load('./steadyState/CaExt001')
                    load('./steadyState/DExt001')
             end
            [~,ind] = min(abs(PARAM.Ca-manyCa));
            Dca = manyD(ind);
            cd(here)
        end
        
        figure(12)
        errDDD = (DDD-Dca)/Dca;
        hold on
        plot(t,errDDD,col)
        grid on
        xlabel('t')
        ylabel('D_{err}')
        title(['Def parameter error Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
    end
    
    %store data for saving
    if plotMovie==1
            Tstore{countShapeLoop-1} = T;
            Ystore{countShapeLoop-1} = Y;
            ElongStore{countShapeLoop-1} = normElon;
            f1Store{countShapeLoop-1} = f1;
            f2Store{countShapeLoop-1} = f2;
            resStore{countShapeLoop-1} = res;
    end

end

%save data
if plotMovie==1
    saveMovieEdgeTracking(Tstore,Ystore,ElongStore,f1Store,f2Store,resStore,color,stepSave,nameNorm,namePhase,nameRes,nameShape,plotLoop,savePlotEnd,solidLine,PARAM,saveDest,saveData);
end

%print simulation time
%display(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











