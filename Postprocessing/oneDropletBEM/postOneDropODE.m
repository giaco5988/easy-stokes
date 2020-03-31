%post processing of test time stepping spectral

close all
clear variables

%data
CaUP = [];
BondUP = 5;
viscUP = 1;
n = 50;
TendUP = 20;
ODEup = 0;          % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 7;             % 1 is extensional flow, 2 is rising droplet
dtUP = 1e-3;
Dup = [0.2 1 0.3];         % initial shape
res = 1;            % results, server or results/forThesis

%option
plotShape = 1;  plotMesh = 0;
plotSnapshot = 0;   timeSnapshot = [0 2.5 5 7.5 9.5];
savePlotShape = 0;  nameSave = 'extensNoBreak';     dest = '/Users/Giacomo/Documents/research_notes/EuroMecBari/presentation/movie/frame/';
getFewShapes = 0;   Tget = linspace(0,18,5);
plotElongNorm = 1;
zoom = 0;
plotRes = 1;
plotPower = 0;
plotAreaNorm = 1;
Tbreak = 2*TendUP;
plotVolErr = 1;
plotCurv = 1;
plotElemNum = 1;
plotPhaseSpacePolar = 0;    dimension = 2;  symmetric = 1;  modes = 4;

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
    name = ['oneDropBEM_ODE=' num2str(ODEup) '_n=' num2str(n) '_BC=' num2str(BC) '_Ca=' num2str(CaUP) '_visc=' num2str(viscUP) '_D=' num2str(Dup) '_maxDT=' num2str(dtUP) '_Tend=' num2str(TendUP) '.mat'];
elseif isempty(CaUP)
    name = ['oneDropBEM_ODE=' num2str(ODEup) '_n=' num2str(n) '_BC=' num2str(BC) '_Bo=' num2str(BondUP) '_visc=' num2str(viscUP) '_D=' num2str(Dup) '_maxDT=' num2str(dtUP) '_Tend=' num2str(TendUP) '.mat'];
end

%upload data
upload = [dir name];
load(upload)

T = allRes{1};
YYY = allRes{2};
VVV = allRes{3};

%number of iteration
ite = round(numel(T));

%initialize
normElon = zeros(ite,1);
elemNum = zeros(ite,1);
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
xGetShape = cell(size(Tget));
yGetShape = cell(size(Tget));
countGet = 1;
indDots = zeros(numel(timeSnapshot),1);

%loops
plotCount = 1;
shift = 0;  islast = 0;
for i = 1:ite
    
   display([num2str(i) ' of ' num2str(ite)])
   
   %get currents modes
   Y = YYY{i};
   x{1} = Y(1:2:end-1)';
   y{1} = Y(2:2:end)';
   xGrid = x{1};
   yGrid = y{1};
   
   elemNum(i) = numel(xGrid);
   
   if plotPhaseSpacePolar==1
            
            %radius modes
            f = LegendreSerie(x,y,modes,symmetric,PARAM);
            f1(i) = f(2);
            f2(i) = f(3);
            f3(i) = f(4);
        
   end
   
   %compute elongation norm
   xcm(i) = centerOfMassBlockAxis(x,y,1,PARAM);
   RforNorm = sqrt((xGrid-xcm(i)).^2+yGrid.^2);
   normElon(i) = abs(RforNorm(end)-xcm(i));
   
   %compute volume abd surface area
   V(i) = axis_int_gauss_vect(xGrid,yGrid);
   Area(i) = surf_gauss_vect(xGrid,yGrid);
   
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
              plot(xGrid,yGrid+shift,'k')
              hold on
              plot(xGrid,-yGrid+shift,'k')
              %grid on
              axis equal
              xlabel('z')
              ylabel('r')
              axis([-5 5 -2 2+3*(totalShapes-1)])

              plotCount = plotCount+1;
              shift = shift-3;
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
   
   %compute curvature
   [nx,ny] = computeNormalVector(x{1},y{1},PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));
   [K1,K2] = computeCurvatureSplines(x{1},y{1},PARAM.orderVariableStokes(1));
   K = K1+K2;
   
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
           
            xGetShape(countGet) = {xGrid};
            yGetShape(countGet) = {yGrid};
            
            countGet = countGet+1;
           
       end
       
   end
   
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
      if BC==7
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
   
   if plotRes==1 || plotPower==1
       
           
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
   
   %stop when there are no more data
   if i~=ite
   if T(i+1)>Tbreak||(T(i+1)==0&&k>1)
       V = V(1:i);
       t = T(1:i);
       res = res(1:i);
       normElon = normElon(1:i);
       disp('Break')
       break;
   end
   end
    
end

if plotElemNum==1

    figure
    plot(t(1:i),elemNum,'ok')
    %xlabel('t')
    %ylabel('\Delta S')
    xyLabelTex('t','N')
    title('Number of elements')
    grid on

end

if plotVolErr==1

    errV = (V-V0)/V0;
    figure
    plot(t,errV)
    xlabel('t')
    ylabel('err_V')
    grid on
    title('Error On Volume')

end

if plotRes==1

    figure
    semilogy(t,res,'k')
    %xlabel('t')
    %ylabel('res')
    xyLabelTex('t','||\mathbf{u} \cdot \mathbf{n}||_\infty')
    title('residuals')
    grid on

end

if plotCurv==1

    figure
    lArc = [0 cumsum(sqrt(diff(xGrid).^2+diff(yGrid).^2))];
    plot(lArc,K1+K2,'k')
    %hold on
    %plot(xGrid,K2)
    %plot(xGrid,K1+K2,'k')
    xyLabelTex('l','k')
    %xlabel('x')
    %ylabel('k')
    %legend('k_1','k_2','k','Location','Best')
    title('curvature')
    grid on

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

if plotElongNorm==1

    figure
    %semilogy(t,normElon,'k')
    plot(t,normElon-1,'k')
    %xlabel('t')
    %ylabel('L/a')
    xyLabelTex('t','L')
    title('Stone elongation')
    grid on
    %axis([0 200 1 10])
    
    if plotSnapshot==1
        
        hold on
        plot(t(indDots),normElon(indDots),'k.','MarkerSize',30)
        
        
    end
    
end

if plotAreaNorm==1

    figure
    semilogy(t(1:i),Anorm(1:i),'k')
    %xlabel('t')
    %ylabel('\Delta S')
    xyLabelTex('t','\Delta S')
    title('Excess area')
    grid on

end

%print simulation time
disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











