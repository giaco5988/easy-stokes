%post processing of a conical motor with a defromable bubble

clear variables
close all

%parameters
res = 0;
BoUP = 0.2;
lambdaUP = 0;
dtUp = .02;
TendUP = 2001;
alphaUP = 1.4;
ODE = 0;
Lup = 20;
addX = 2;
%nDropUP = 80;   nWallUP = 10;
nDropUP = round(alphaUP*50)+10;   nWallUP = 10;
nElemUP = nDropUP+nWallUP*Lup+2*nWallUP;
x0Upload = 0;
repUP = 5;

color = get(gca,'ColorOrder');

%options
plotMotorUp = 1;    plotIteUp = 1;     savePlot = 0;    plotMesh = 0;
threeD = 0; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.3];
plotRes = 0;
plotDropVel = 1;
plotFilm = 1;
plotBubbleArea = 0;
plotBubbleVolume = 1;
plotBubbleXcm = 0;
plotLastCurvature = 1;
plotVelField = 0;   Tplot = 300;    substitutePoint = 1;    coeffSub = 1; noIn = 0;
plotVelFieldSlice = 0;  xSlice = 2; nSlice = 100;   Ytop = 1.4;
cutX = 5;   shift = 5;  cutY = 0.99;   nx = 41;    ny = 40;    postDropFrame = 0;
plotVelDecay = 0;   nDecay = 5;

%filename
filename = ['verticalTube_ODE=' num2str(ODE) '_rep=' num2str(repUP) '_x0=' num2str(x0Upload) '_alpha=' num2str(alphaUP) '_Bond=' num2str(BoUP) '_lambda=' num2str(lambdaUP) '_L=' num2str(Lup) '_nStart=' num2str(nElemUP) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];

%source
if res==1
   source = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==0
   source = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==2
   source = '~/Documents/MATLAB/test/verticalPipeDrop/resTest/';
end

BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/verticalTube';
here = pwd;

saveDest = '/Users/Giacomo/Documents/research_notes/verticalPipe/movies/frames/';
namePlot = 'plotVerticalPipe';

%load file
load([source filename])
PARAM.kernelFreeSpace = 1;
%PARAM.ellipseShape = 1;
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
tParametricPost = {linspace(0,1,nWallUP) linspace(0,1,nWallUP*Lup) linspace(0,1,nWallUP) nan};
for i = 1:loopUp
    
   display([num2str(i) ' of ' num2str(loopUp)])
    
   %current shape and motor position
   Yhere = Y{i};
   if ODE==0
    Vhere = V{i};
   end
   xBubble = Yhere(1:2:end-1);
   yBubble = Yhere(2:2:end);
   
   %current location   TEMPORARY!!!
   PARAM_now = PARAM;
   
   %compute curent shape
   [x,y] = buildGeometryPanelsParametricPost(tParametricPost,PARAM_now);
   x{4} = xBubble';
   y{4} = yBubble';
   
   %residuals
   if plotRes==1 || plotDropVel==1
    if ODE~=0
        cd(BEM)
        Vhere = computeVelocityRising(T(i),Yhere,tParametricPost,PARAM_now);
        cd(here)
    end
    [nxBubble,nyBubble] = computeNormalVector(xBubble',yBubble',PARAM.orderVariableStokes(4),PARAM.orderGeometryStokes(4),PARAM.SPlinesType(4));
    dropVel(i) = DropVelocityAxis(xBubble',yBubble',Vhere(1:2:end-1).*nxBubble'+Vhere(2:2:end).*nyBubble');
    frontVel(i) = Vhere(1);
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
   if plotMotorUp==1 && sum(i==(1:plotIteUp:numel(T)))
        
        %block coordinates
        figure(1)
        plotGeometryStokes(x,y,threeD,theta3D,color3D,transp,0,PARAM_now)
        %plot(xBubble(1),yBubble(1),'or')
        xcmPlot = center_mass(xBubble,yBubble);
        axis([xcmPlot-4 xcmPlot+4 -1 1])
        grid on
        hold off
        xlabel('z','Interpreter','latex')
        ylabel('r','Interpreter','latex')
        %axis off
        title(['Bo=' num2str(BoUP)],'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        
        if plotMesh==1
            
            hold on
            plot(x{4},y{4},'rx')
            hold off
            
        end
        drawnow
        
        %subplot(2,1,2)
%         figure(1)
%         plot(T(1:i),Ucone(1:i))
%         hold on
%         plot(T(i),Ucone(i),'.','MarkerSize',40,'Color',color(1,:))
%         grid on
%         xlabel('t')
%         ylabel('U_{cone}')
%         hold off
%         axis([0 500 -2.5*1e-3 -0.5*1e-3])
%         %title('Cone Velocity')
%         drawnow
        
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
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel,Vdrop] = computeVelocityRising(T(i),Y{i},tParametricBase,PARAM);
       cd(here)
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
   
   %plot velocity field
   if plotVelFieldSlice==1 && T(i)==Tplot
       
       %create grid
       Ygrid = linspace(0,Ytop,nSlice);
       Xgrid = xSlice*ones(1,numel(Ygrid));
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,PARAM);
       cd(here)
       
       %on the right of the motor
       [~,~,uxSlice,uySlice] = computeVelPressField(Xgrid,Ygrid,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub);
       
       figure
       subplot(2,1,1)
       plotGeometryStokes(x,y,threeD,theta3D,color3D,transp,PARAM_now)
       axis([-5 15 -2 2])
       grid on
       title(['t=' num2str(T(i))])
       drawnow
       plot(Xgrid,Ygrid,'r')
       
       subplot(2,1,2)
       plot(Ygrid,uxSlice)
       hold on
       grid on
       plot(Ygrid,uySlice)
       xlabel('y')
       ylabel('u')
       
       %compute flow rate
       Qflow = 2*pi*trapz(Ygrid,Ygrid.*uxSlice);
       display(Qflow)
       
   end
   
   %plot velocity field
   if plotVelDecay==1 && T(i)==Tplot
       
       %create grid
       Xgrid = linspace(-cutX,cutX,nx)+shift;
       Ygrid = linspace(0,cutY,ny);
       [Xgrid,Ygrid] = meshgrid(Xgrid,Ygrid);
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,PARAM);
       cd(here)
       
       %on the right of the motor
       Xright = logspace(1,nDecay,100);
       Yright = zeros(1,numel(Xright));
       [~,~,uxRight,uyRight] = computeVelPressField(Xright,Yright,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub);
       
       %on the left of the motor
       Xleft = -logspace(1,nDecay,100);
       Yleft = zeros(1,numel(Xleft));
       [~,~,uxLeft,uyLeft] = computeVelPressField(Xleft,Yleft,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub);
       
       %upward
       Yupward = logspace(1,nDecay,100);
       Xupward = zeros(1,numel(Yupward));
       [~,~,uxUpward,uyUpward] = computeVelPressField(Xupward,Yupward,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub);
       
       figure
       loglog(Xright,sqrt(uxRight.^2+uyRight.^2),'-x')
       hold on
       loglog(Yupward,sqrt(uxUpward.^2+uyUpward.^2),'-x')
       loglog(-Xleft,sqrt(uxLeft.^2+uyLeft.^2),'-x')
       grid on
       xlabel('\rho')
       ylabel('|U|')
       title('Velocity decay')
       drawnow
       
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
    %hold on
    plot(T(1:i),hFilmPost(1:i))
    hold on
    %plot(T(1:i),minDL(1:i))
    grid on
    xyLabelTex('t','h')
    title('film thickness')
    %legend('Film thickness','Minimum element size','Location','Best')
    drawnow
    
    figure(10)
    %hold on
    plot(xcm(1:i),hFilmPost(1:i))
    grid on
    xyLabelTex('z_{cm}','h')
    title('film thickness')
    %legend('Film thickness','Minimum element size','Location','Best')
    drawnow
       
    figure(10)
    %hold on
    plot(maxBubble(1:i),maxBubble(1:i)-minBubble(1:i),'k')
    grid on
    xyLabelTex('z_{front}','\Delta l')
    title('Bubble lenght')
    %legend('Film thickness','Minimum element size','Location','Best')
    drawnow
    
end

%plot motor velocity
if plotRes==1
    
    figure(8)
    %hold on
    semilogy(T(1:i),manyRes(1:i))
    %plot(T(1:i),manyRes(1:i))
    grid on
    xyLabelTex('t','R')
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














