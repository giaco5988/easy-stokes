%post processing of test time stepping spectral

close all
clear variables

%set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%data
manyVisc = 1*ones(26,1);
nnn = 50*ones(5,1);
%manyCa = [0.125 0.0625 0.0625 0.0625 0.0625 0.0625 0.0625];
manyCa = [-0.4 -.35*ones(2,1)' -.38*ones(2,1)'];
manyTend = [100 100+25 100+26 100+12 100+13];
%manyTend = [70 100 251.728 251.73 84];
%manyTend = 251.7:2e-3:251.75;
ODE = 2;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1;         % 1 is extensional flow, 2 is rising droplet
volCorr = 0;
DT = 1e-2*ones(26,1);
manyLegendre = 1*ones(26,1);
manyD = 0*ones(26,1);          % initial shape
CPUs = 16;
res = 1;        % results or server

%plot shape at these times
plotSnapshots = 1;
% cellShapes{1} = [];
% cellShapes{2} = [51 64 80 98];
% cellShapes{3} = [51 75 100 125 149]+0.728;
% cellShapes{4} = [51 72 95 119]+0.73;
% cellShapes{5} = [52 60 70 81.75];
cellShapes{1} = [];
cellShapes{2} = [25 50 75 100];
cellShapes{3} = [26 50 75 100];
cellShapes{4} = [12 50 75 100];
cellShapes{5} = [13 50 75 100];


%legend
cellLegend = cell(size(nnn));

%option
plotShape = 1;
plotElong = 1;
plotExcessArea = 0;
plotVolumeErr = 0;
plotPhaseSpacePolar = 1;    dimension = 2;  symmetric = 0;  modes = 7;  beforeEnd = 8;
plotPhaseSpacePOF = 0;

%saving options
plotMovie = 0;
saveData = [0 0 0];
plotLoop = 2;
stepSave = [10 2 10 10 2];
savePlotEnd = [100 100 150 120 90];
solidLine = [1 1 1 1 1];
nameShape = 'shapeStoneTogether';
nameNorm = 'normStoneTogether';
namePhase = 'phaseStoneTogether';
nameSI = 'SIprf';

%destination for saving plot
saveDest = '/Users/Giacomo/Documents/research_notes/edge_state/draftPRF/movie/frame/';

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

%get colot rgb
color = get(gca,'ColorOrder');
close gcf

%iintialize
if plotMovie==1
    
   Tstore = cell(size(nnn));
   Ystore = cell(size(nnn));
   ElongStore = cell(size(nnn));
   f1Store = cell(size(nnn));
   f2Store = cell(size(nnn));
   resStore = cell(size(nnn));
    
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
    Ca = manyCa(l);
    legendre = manyLegendre(l);
    if plotSnapshots==1
        Tshape = cellShapes{l};
    end
    
        %store element in cell for legend
        cellLegend(l) = {['Shape ' num2str(l)]};

        %filename
        name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(D) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '.mat'];

        %upload data
        upload = [dir name];
        load(upload)

        %number of iteration
        ite = numel(T);

        %initialize chebfun
        if PARAM.legendre==0
            x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
        end

        %initialize
        res = zeros(ite,1);
        V = zeros(ite,1);
        elongation = zeros(ite,1);
        excessArea = zeros(ite,1);
        V0 = PARAM.V0;
        f1 = zeros(ite,1);
        f2 = zeros(ite,1);
        f3 = zeros(ite,1);
        if plotSnapshots==1
            indDots = zeros(numel(Tshape),1);
        end

        %loops
        countShift = 1;
        countShape = 1;
        for i = 1:ite

           %display([num2str(i) ' of ' num2str(ite)])

           %get currents modes
           xyMode = Y(i,:)';
           xMode = Y(i,1:2:end-1)';
           yMode = Y(i,2:2:end)';
           

           %from modes to grid
           [x,y] = fromModesToGrid(xMode,yMode,PARAM);
           
           %defromation parameter
%            L = max(x)-min(x);
%            B = 2*y((PARAM.n+1)/2);
%            Dellipse(i) = (L-B)/(L+B);
%            delta(i) = L/B-1;

           if plotPhaseSpacePolar==1
            
                %radius modes
                f = LegendreSerie(x,y,modes,symmetric,PARAM);
                f1(i) = f(3);
                f2(i) = f(5);
                f3(i) = f(7);
%                 f1(i) = yMode(3);
%                 f2(i) = yMode(5);
%                 f3(i) = f(7);
        
           end
    
           if 1 %Ca>=0
                elongation(i) = x(1);
           else
                [~,ind] = max(y);
                elongation(i) = y(ind);
           end

           %compute volume
           V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);

           R0 = nthroot(V(i)/4/pi*3,3);
           excessArea(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM) - 4*pi*R0^2;

           t = T;

           if plotShape==1

              xMid = (max(x)+min(x))/2;

              figure(1) 
              if BC==2
                    plot(y,-x,'k-')
                    hold on
                    grid on
                    plot(-y,-x,'k')
                    title('Rising Droplet')
                    axis([-2 2 -2 2])
              elseif BC==1
                    plot(x,y,'k-')
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
           
           if plotSnapshots==1 && sum(abs(T(i)-Tshape)<1e-4)
              
              if Ca>=0
                shiftX = 10;
                shiftY = 2;
              else
                  shiftX = 4;
                  shiftY = 5;
              end
              
              figure(12)
              plot([x flip(x)]+shiftX*(l-1),[y -flip(y)]+shiftY*(countShift-1),'Color',color(l,:))
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

           %stop when there are no more data
           if i~=ite
               if T(i+1)>Tbreak||(T(i+1)==0&&i>1)
                   ForceFree = ForceFree(1:i);
                   V = V(1:i);
                   t = T(1:i);
                   res = res(1:i);
                   excessArea = excessArea(1:i);
                   elongation = elongation(1:i);
                   disp('Break')
                   break;
               end
           end

        end
        
        if plotElong==1

            %plot deformation parameter semilogy
            figure(2)
            if l>1
                hold on
            end
            if Ca>=0
                semilogy(t(1:end-beforeEnd),elongation(1:end-beforeEnd))
            else
                semilogy(t(1:end-beforeEnd),elongation(1:end-beforeEnd))
            end
            %xlabel('t')
            %ylabel('max(x)')
            xyLabelTex('t','L')
            grid on
            if Ca>=0
                axis([0 150 1 10])
            else
                %axis([0 150 1 10])
            end
            title('Elongation')
            
            if plotSnapshots
                
               hold on
               semilogy(t(indDots),elongation(indDots),'.k','MarkerSize',30) 
                
            end
        
        end
        
        if plotPhaseSpacePolar==1
        
            figure(11)

            if dimension==2

                hold on
                plot(f1(1:end-beforeEnd),f2(1:end-beforeEnd))
                grid on

                if plotSnapshots
                
                    hold on
                    plot(f1(indDots),f2(indDots),'.k','MarkerSize',30) 
                
                end

            elseif dimension==3

                hold on
                plot3(f1(1:end-beforeEnd),f2(1:end-beforeEnd),f3(1:end-beforeEnd))
                zlabel('f_3')
                grid on
                
                if plotSnapshots
                
                    hold on
                    plot3(f1(indDots),f2(indDots),f3(indDots),'.k','MarkerSize',30) 
                
                end

            end

            xlabel('f_2')
            ylabel('f_4')
            title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
        end
        
        if plotExcessArea==1

            %plot excess surface area semilogy
            figure(3)
            if l>1
                hold on
            end
            semilogy(t,excessArea)
            xlabel('t')
            ylabel('\Delta A')
            grid on
            %axis([0 8 0 1])
            title('Excess area')
        
        end
        
        if plotPhaseSpacePOF==1
            
            %compute lDot
            lDot = diff(elongation(1:end-beforeEnd))./diff(t(1:end-beforeEnd));
            lDot = [lDot; lDot(end)];
        
            figure(5)
            hold on
            plot(elongation(1:end-beforeEnd),lDot)
            grid on

            if plotSnapshots
                
                hold on
                plot(elongation(indDots),lDot(indDots),'.k','MarkerSize',30) 
                
            end

            xlabel('l')
            ylabel('dl/dt')
            title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
        end

        if plotVolumeErr
            figure(4)
            if l>1
                hold on
            end
            errV = (V-V0)/V0;
            plot(t,errV)
            xlabel('t/\tau_r')
            ylabel('err_V')
            grid on
            title('Error on volume')
        end
        
        %print simulation time
        display(['Simulation time T=' num2str(simulationTime/60) ' minutes'])
        
        %store data for saving
        if plotMovie==1
            Tstore{l} = T;
            Ystore{l} = Y;
            ElongStore{l} = elongation;
            f1Store{l} = f1;
            f2Store{l} = f2;
            %resStore{l} = res;
        end
    
end

%save data
if plotMovie==1
    %saveMovieCompareExtensional(Tstore,Ystore,ElongStore,f1Store,f2Store,color,stepSave,nameNorm,namePhase,nameShape,plotLoop,savePlotEnd,solidLine,PARAM,saveDest,saveData);
    %saveMovieCompareExtensionalHardCoded(Tstore,Ystore,ElongStore,f1Store,f2Store,color,stepSave,nameNorm,namePhase,nameShape,plotLoop,savePlotEnd,solidLine,PARAM,saveDest,saveData);
    saveMovieCompareExtensionalHardCodedSameTime(Tstore,Ystore,ElongStore,f1Store,f2Store,color,stepSave,nameSI,plotLoop,savePlotEnd,solidLine,PARAM,saveDest,0);
end

if plotElong
    figure(2)
    legend(cellLegend,'Location','Best')
end

if plotExcessArea
    figure(3)
    legend(cellLegend,'Location','Best')
end

if plotVolumeErr
    figure(4)
    legend(cellLegend,'Location','Best')
end






