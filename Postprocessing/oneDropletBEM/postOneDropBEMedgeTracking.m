%post processing of test time stepping spectral

close all
clear variables

%data
CaUP = [];
BondUP = 1.5;
viscUP = 1;
nUP = 100;
%TendUP = 20;
ODEup = 0;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BCup = 7;         % 1 is extensional flow, 2 is rising droplet
dtUP = 1e-3;
edgeLoopUP = 1002;
deltaUP = 1e-6;
spanPostLoop = 1:edgeLoopUP;

res = 0;        %results, server or results/forThesis

%option
plotShape = 0;  plotMesh = 0;
plotSnapshot = 0;   timeSnapshot = [0 10 20 26];
savePlotShape = 0;  nameSave = 'extensNoBreak';     dest = '/Users/Giacomo/Documents/research_notes/EuroMecBari/presentation/movie/frame/';
plotElongNorm = 1;
plotRes = 1;
plotAreaNorm = 0;
plotVolErr = 0;
plotPhaseSpacePolar = 1;    dimension = 2;  symmetric = 0;  modes = 4;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==2
    dir = '~/Documents/MATLAB/droplet_simulations/results/forThesis/manyExtensionalFlow/';
end

%filename
if isempty(BondUP)
    %name = ['oneDropBEM_ODE=' num2str(ODEup) '_n=' num2str(n) '_BC=' num2str(BC) '_Ca=' num2str(CaUP) '_visc=' num2str(viscUP) '_D=' num2str(Dup) '_maxDT=' num2str(dtUP) '_Tend=' num2str(TendUP) '.mat'];
elseif isempty(CaUP)
    name = ['edgeTrackingDropBEM_edgeLoop=' num2str(edgeLoopUP) '_deltaEdge=' num2str(deltaUP) '_ODE=' num2str(ODEup) '_BC=' num2str(BCup) '_Bo=' num2str(BondUP) '_visc=' num2str(viscUP) '_n=' num2str(nUP) '_maxDT=' num2str(dtUP) '.mat'];
end

%upload data
upload = [dir name];
load(upload)

%loop on edge tracking
for k = spanPostLoop
    
display([num2str(k) ' of ' num2str(edgeLoopUP)])

T = Tedge{k};
YYY = Yedge{k};
VVV = Vedge{k};

if isempty(T)
        %disp('No data')
        break;
end

%number of iteration
ite = round(numel(T));

%initialize
normElon = zeros(ite,1);
res = zeros(ite,1);
V = zeros(ite,1);
Area = zeros(ite,1);
Evisc = zeros(ite,1);
Dellipse = zeros(ite,1);
xcm = zeros(ite,1);
Anorm = zeros(ite,1);
f1 = zeros(ite,1);
f2 = zeros(ite,1);
f3 = zeros(ite,1);
V0 = 4/3*pi;
countGet = 1;
indDots = zeros(numel(timeSnapshot),1);

%loops
plotCount = 1;
shift = 0;  islast = 0;
for i = 1:ite
    
   %get currents modes
   Y = YYY{i};
   x{1} = Y(1:2:end-1)';
   y{1} = Y(2:2:end)';
   xGrid = x{1};
   yGrid = y{1};
   
   if plotPhaseSpacePolar==1
            
            %radius modes
            f = LegendreSerie(x{1}',y{1}',modes,symmetric,PARAM);
            f1(i) = f(2);
            f2(i) = f(3);
            f3(i) = f(4);
        
   end
   
   %compute elongation norm
   RforNorm = sqrt(xGrid.^2+yGrid.^2);
   %xcm(i) = centerOfMassBlockAxis(x,y,1,PARAM);
   xcm(i) = center_mass(x{1},y{1});
   %normElon(i) = abs(RforNorm(end)-xcm(i));
   normElon(i) = abs(xGrid(end)-xcm(i));
   
   %compute volume abd surface area
   V(i) = axis_int_gauss_vect(xGrid,yGrid);
   Area(i) = surf_gauss_vect(xGrid,yGrid);
   
   if plotAreaNorm==1
       
       R0 = nthroot(V(i)/4/pi*3,3);
       %Anorm(i) = (Area(i) - 4*pi*R0^2)./(Area(1) - 4*pi*R0^2);
       Anorm(i) = Area(i) - 4*pi*R0^2;
       
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
   
   %defromation parameter
   L = max(xGrid)-min(xGrid);
   B = 2*yGrid(round((PARAM.n+1)/2));
   Dellipse(i) = (L-B)/(L+B);
   
   if PARAM.typeBCstokes==7

        OutIn  = abs(Dellipse(i)) > PARAM.convergeShapeEdge;

   end
   
   if OutIn==1
        col = '-';
    elseif OutIn==0
        col = '--';
    end
   
   %get time
   t = T;
   
   if plotShape==1
       
      xMid = (max(xGrid)+min(xGrid))/2;
       
      figure(1) 
      plotGeometryStokes(x,y,0,[],[],[],0,PARAM)
      %xlabel('r')
      %ylabel('x')
      xyLabelTex('z','r')
      hold off
      title('Drop shape')
      
      if plotMesh==1
         
          hold on
          plot(x{1},y{1},'xr')
          hold off
          
      end
      axis equal
      axis([-4 2 -2 2])
      if BCup==7
          camroll(90)
      end
      drawnow
      
      if savePlotShape==1
          
        disp('Save')
        grid off
        axis off
        title('')
        print('-dpng','-loose','-r400',[dest nameSave sprintf('%03d',i-2) '.png'])
          
      end
      
      if zoom==1
          
         axis([-0.1 0.1 -1.2 -1.0]) 
          
      end
      
      drawnow
      
   end
   
   if plotRes==1
       
           
           if PARAM.ODE==0
               
               Vhere = VVV{i};
               UxN = Vhere(1:2:end-1);
               UyN = Vhere(2:2:end);
               UnABS = sqrt(UxN.^2+UyN.^2);
               res(i) = norm(UnABS,Inf);
               
           else
               error('Not implemented')
           end
           
%            [u,v] = fromModesToGrid(uvModes(1:2:end-1),uvModes(2:2:end),PARAM);
%            fE = fx.*u + fy.*v;
%            Evisc(i) = axisIntSpectralCurvilinear(x,y,fE,PARAM);
       
   end
   
%    %stop when there are no more data
%    if i~=ite
%    if T(i+1)>Tbreak||(T(i+1)==0&&k>1)
%        V = V(1:i);
%        t = T(1:i);
%        res = res(1:i);
%        normElon = normElon(1:i);
%        disp('Break')
%        break;
%    end
%    end
    
end

if plotVolErr==1

    errV = (V-V0)/V0;
    figure(2)
    if k>1
        hold on
    end
    plot(t,errV)
    xlabel('t')
    ylabel('err_V')
    grid on
    title('Error On Volume')

end

if plotRes==1

    figure(3)
    semilogy(t,res,col)
    hold on
    %xlabel('t')
    %ylabel('res')
    xyLabelTex('t','||\mathbf{u} \cdot \mathbf{n}||_\infty')
    title('residuals')
    grid on

end

%plot phase space
if plotPhaseSpacePolar==1
        
        figure(11)
        
        if dimension==2
            
            hold on
            plot(f1(1:i),f2(1:i),col)
            grid on
            
            if plotSnapshot==1
                plot(f1(indDots),f2(indDots),'k.','MarkerSize',30)
            end
            
        elseif dimension==3
            
            hold on
            plot3(f1(1:i),f2(1:i),f3(1:i),col)
            zlabel('f_3')
            grid on
            
        end
        
        xlabel('f_1')
        ylabel('f_2')
        %title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        title(['Phase space Bo=' num2str(BondUP) ' \lambda=' num2str(PARAM.visc)])
        
end

if plotElongNorm==1

    figure(4)
    %semilogy(t,normElon,'k')
    if k>1
        hold on
    end
    plot(t,normElon-1,col)
    %semilogy(t,normElon-1,col)
    %xlabel('t')
    %ylabel('L/a')
    xyLabelTex('t','L')
    title('Elongation')
    grid on
    %axis([0 200 1 10])
    
    if plotSnapshot==1
        
        hold on
        plot(t(indDots),normElon(indDots),'k.','MarkerSize',30)
        
        
    end
    
end

if plotAreaNorm==1

    figure(5)
    if k>1
        hold on
    end
    plot(t(1:i),Anorm(1:i),col)
    %xlabel('t')
    %ylabel('\Delta S')
    xyLabelTex('t','\Delta S')
    title('Excess area')
    grid on

end

end

%print simulation time
try
    disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])
catch
    disp('Simulation has not yet finished')
end











