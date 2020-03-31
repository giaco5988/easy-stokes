%post processing of a conical motor with a defromable bubble

clear variables
close all

color = get(gca,'ColorOrder');

%% Add libraries paths
REPOSITORY_NAME = '~/Documents/MATLAB/';    % path to the repository
add_paths_tutorials(REPOSITORY_NAME);

%% Uplaod parameters
BoUP = 0;           % Bond number
CaUP = 0.05;        % Capillary number
lambdaUP = 10;      % viscosity ratio
dtUp = .1;          % time step
TendUP = 500;       % total simulated time
alphaUP = 1.1;      % drop size
Lup = 10;           % length of the capillary
nElemUP = 160;      % total number of elements (at start)
x0Upload = 0;       % initial drop position
repUP = 0;          % repulsive forces option

%% Plotting options
plotBubble = 1;    plotIteUp = 1;     savePlot = 0;    plotMesh = 0;                    % if 1, do plot and choose options
threeD = 0; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.3]; % make plot 3D
plotRes = 1;                                                                            % plot residuals
plotDropVel = 0;                                                                        % plot drop velocity
plotFilm = 0;                                                                           % plot film thickness
plotBubbleArea = 0;                                                                     % plot bubble area
plotBubbleVolume = 0;                                                                   % plot bubble volume
plotBubbleXcm = 0;                                                                      % plot center of mass position
plotLastCurvature = 0;                                                                  % plot curcature
plotVelField = 0;   Tplot = 450;    substitutePoint = 1;    coeffSub = 0.5; noIn = 0;   % plot velocity field
cutX = 5;   shift = 0;  cutY = 0.99;   nx = 100;    ny = 10;    postDropFrame = 1;

%filename
filename = ['verticalTube_ODE=' num2str(0) '_rep=' num2str(repUP) '_x0=' num2str(x0Upload) '_alpha=' num2str(alphaUP) '_Bond=' num2str(BoUP) '_Ca=' num2str(CaUP) '_lambda=' num2str(lambdaUP) '_L=' num2str(Lup) '_nStart=' num2str(nElemUP) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];

%source
source = '../tutorial_results/';

%load file
load([source filename])
PARAM.kernelFreeSpace = 1;
PARAM.ellipseShape = [0 0 0 1];
PARAM.D = [0 0 0 0];

%number
if T(end)==0
    loopUp = find(T==0,2);
    loopUp = loopUp(2)-1;
else
    loopUp = numel(T);
end
%initialize
manyRes = zeros(numel(T),1);
dropVel = zeros(numel(T),1);
frontVel = zeros(numel(T),1);
Area = zeros(numel(T),1);
Volume = zeros(numel(T),1);
xcm = zeros(numel(T),1);
maxBubble = zeros(numel(T),1);
minBubble = zeros(numel(T),1);
hFilmPost = zeros(numel(T),1);
minDL = zeros(numel(T),1);
tParametricPost = {linspace(0,1,PARAM.n(1)+1) linspace(0,1,PARAM.n(2)+1) linspace(0,1,PARAM.n(3)+1) nan};
for i = 1:loopUp
    
   display([num2str(i) ' of ' num2str(loopUp)])
    
   %current shape and motor position
   Yhere = Y{i};
   Vhere = V{i};
   xBubble = Yhere(1:2:end-1);
   yBubble = Yhere(2:2:end);
   
   %current location
   PARAM_now = PARAM;
   
   %compute curent shape
   [x,y] = buildGeometryPanelsParametricPost(tParametricPost,PARAM_now);
   x{4} = xBubble';
   y{4} = yBubble';
   
   %residuals
   if plotRes==1 || plotDropVel==1
    [Vhere,~,~,~,~,dropVel(i)] = computeVelocityRising(T(i),Yhere,tParametricPost,PARAM_now);
    [nxBubble,nyBubble] = computeNormalVector(xBubble',yBubble',PARAM.orderVariableStokes(4),PARAM.orderGeometryStokes(4),PARAM.SPlinesType(4));
%     dropVel(i) = DropVelocityAxis(xBubble',yBubble',Vhere(1:2:end-1).*nxBubble'+Vhere(2:2:end).*nyBubble');
    frontVel(i) = Vhere(1);
    if PARAM.dropFrame~=0
        frontVel(i) = frontVel(i) + dropVel(i);
    end
    res = NormalVelocityDropFrame(Yhere(1:2:end-1)',Yhere(2:2:end)',Vhere(1:2:end-1),Vhere(2:2:end));
    manyRes(i) = norm(res,Inf);
   end
   
   %bubble surface area
   Area(i) = surf_gauss_vect(xBubble',yBubble');
   
   %bubble volume
   Volume(i) = axis_int_gauss_vect(xBubble',yBubble');
   
   %bubble center of mass
   xcm(i) = center_mass(xBubble',yBubble');
   maxBubble(i) = max(xBubble);
   minBubble(i) = min(xBubble);
   %xcm(i) = centerOfMassBlockAxis(x,y,2,PARAM);
   
   %compute distance wall-drop
   %hFilmPost(i) = min(distWallDrop(x{4},y{4},x{2},y{2}));
   hFilmPost(i) = min(1-y{4});
   dx = diff(xBubble);
   dy = diff(yBubble);
   dl = sqrt(dx.^2+dy.^2);
   minDL(i) = min(dl);
   
   %plot shape
   if plotBubble==1 && sum(i==(1:plotIteUp:numel(T)))
        
        %block coordinates
        figure(1)
        plotGeometryStokes(x,y,threeD,theta3D,color3D,transp,0,PARAM_now)
        %plot(xBubble(1),yBubble(1),'or')
        xcmPlot = center_mass(xBubble,yBubble);
        plot(xcmPlot, 0, 'or')
        grid on
        hold off
        xlabel('z','Interpreter','latex')
        ylabel('r','Interpreter','latex')
        %axis off
        title(['Bo=' num2str(BoUP) ' Ca=' num2str(CaUP)],'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        
        if plotMesh==1
            
            hold on
            plot(x{4},y{4},'rx')
            hold off
            
        end
        drawnow
        
        if savePlot==1
            grid off
            print('-dpng','-loose','-r100',[saveDest namePlot sprintf('%03d',round((i-1)/plotIteUp)) '.png'])
        end
        
   end
   
   %plot velocity field
   if plotVelField==1 && T(i)==Tplot
       
       %create grid
       Xgrid = linspace(-cutX,cutX,nx)+shift;
       Ygrid = linspace(0,cutY,ny);
       [Xgrid,Ygrid] = meshgrid(Xgrid,Ygrid);
       
       %compute
       [~,yStokes,xVel,yVel,PARAMvel,Vdrop] = computeVelocityRising(T(i),Y{i},tParametricBase,PARAM);
       [Xsing,Ysing,ux,uy] = computeVelPressField(Xgrid,Ygrid,xVel,yVel,yStokes,[0 0],0,PARAMvel,substitutePoint,coeffSub,noIn);
       frame = 'lab';
       if postDropFrame==1
        ux = ux-Vdrop;
        frame = 'drop';
       end
       
       figure
       [~,h1] = contourf(Xsing,Ysing,sqrt(ux.^2+uy.^2),500);
       hold on
       [~,h2] = contourf(Xsing,-Ysing,sqrt(ux.^2+uy.^2),500);
       hStream1 = streamslice(Xsing,Ysing,ux,uy,0.2);
       hStream2 = streamslice(Xsing,-Ysing,ux,-uy,0.2);
       xlabel('z')
       ylabel('r')
       colorbar
       title(['Velocity field in ' frame ' frame'])

       set(hStream1,'LineWidth',1,'Color','white')
       set(hStream2,'LineWidth',1,'Color','white')
       set(h1,'LineColor','none')
       set(h2,'LineColor','none')
       hold on
       plotGeometryStokes(x,y,threeD,theta3D,color3D,transp,0,PARAM_now)
       axis([-cutX cutX -cutY cutY])
       
       
   end
   
end

%derivative
D1 = finiteDifference1D(i,[2 0],1);

%area for a sphere
Rarea = nthroot(3/4/pi*Volume,3);
AreaSphere = 4*pi*Rarea;

if plotLastCurvature==1
    
   %curvature
   PARAM_now.n(4) = numel(x{4})-1;
   %[nx,ny] = computeNormalVector(x{4},y{4},PARAM_now.orderVariableStokes(4),PARAM_now.orderGeometryStokes(4),PARAM_now.SPlinesType(4));
   [k1,k2] = computeCurvatureSplines(x{4},y{4},PARAM_now.orderVariableStokes(4));
   k = k1+k2;
    
   figure(20)
   dx = diff(x{4}); dy = diff(y{4}); dl = sqrt(dx.^2+dy.^2);
   %lArc = [0 cumsum(dl)];
   %plot(lArc,k,'-x')
   plot(x{4},k,'-k')
   xyLabelTex('z','K')
   title('Curvature')
   grid on
    
end

%plot bubble volume and its variation
if plotBubbleVolume==1
    
    figure(3)
    %subplot(2,1,1)
    %hold on
    Vstart = 4/3*pi*alphaUP^3;
    err = abs(Volume(1:i)-Vstart)/Vstart;
    semilogy(T(1:i),err)
    grid on
    xyLabelTex('t','errV')
    title('Error on Volume')
    drawnow
      
end

%plot bubble surface area
if plotBubbleArea==1
    
    figure
    %subplot(2,1,1)
    plot(T(1:i),Area(1:i)-AreaSphere(1:i),'-k')
    grid on
    xyLabelTex('t','A')
    title('Excess surface area')
    drawnow
    
    dAdt = (D1*Area(1:i))./(D1*T(1:i));
    dAdtSphere = (D1*AreaSphere(1:i))./(D1*T(1:i));
    
    %subplot(2,1,2)
    figure
    plot(T(1:i),dAdt(1:i)-dAdtSphere(1:i),'-k')
    grid on
    xyLabelTex('t','dA/dt')
    title('Surface area variation')
    drawnow
      
end

%plot bubble center of mass
if plotBubbleXcm==1
    
    figure
    %subplot(2,1,1)
    plot(T(1:i),xcm(1:i),'-k')
    grid on
    xyLabelTex('t','z_{cm}')
    title('Bubble position')
    drawnow
      
end

%plot film thicjness
if plotFilm==1
    
    figure(7)
    plot(T(1:i),hFilmPost(1:i))
    hold on
    grid on
    xyLabelTex('t','h')
    title('film thickness')
    drawnow
    
end

%plot motor velocity
if plotRes==1
    
    figure(8)
    semilogy(T(1:i),manyRes(1:i))
    grid on
    xyLabelTex('t','|\bf u \cdot n|')
    title('Residuals')
    drawnow
       
end

if plotDropVel==1
    
    figure(9)
    plot(T(1:i),dropVel(1:i),'k')
    grid on
    xyLabelTex('t','u_d')
    title('Droplet velocity')
    drawnow
    
    figure
    %subplot(2,1,1)
    plot(maxBubble(1:i),dropVel(1:i),'-k')
    grid on
    xyLabelTex('z_{front}','u_d')
    title('Bubble velocity')
    drawnow
    
    figure
    %subplot(2,1,1)
    plot(maxBubble(1:i),frontVel(1:i),'-k')
    grid on
    xyLabelTex('z_{front}','u_{front}')
    title('Front velocity')
    drawnow
       
end

try
    disp(['Simulation time is ' num2str(simulationTime/60/60) ' hours'])
catch
    disp('Simulation has not yet finished')
end














