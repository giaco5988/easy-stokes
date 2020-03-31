%post processing of a conical motor with a defromable bubble

clear variables
close all

%parameters
res = 0;
CaUP = 1e-2;
%thetaUp = linspace(0,pi/32,17);
%thetaUp = thetaUp(10);
thetaUp = 1;
dtUp = 2e-2;
nElemUP = 124;
TendUP = 5000;
coeffRepUP = -10;
alphaUP = 0.5;
ODE = 0;
repTypeUP = 8;
massFlow = 1;
x0upload = 5;

color = get(gca,'ColorOrder');

%options
plotMotorUp = 0;    plotIteUp = 5;     savePlot = 0;
threeD = 0; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.3];
plotConeVel = 1;
stopSphere = 1; stopValue = 2e-3;
plotBubbleArea = 0;
plotBubbleVolume = 0;
plotBubbleXcm = 0;
plotVelField = 1;   Tplot = [0 2000 2500 2600];    eliminateBlockPostprocessing = [0 0];    streamslicePLot = 0;
substitutePoint = 0;    coeffSub = 0.5;     noIN = 0;
plotVelFieldSlice = 0;  xSlice = 2; nSlice = 100;   Ytop = 1.4;
cutX = 8;   shift = 5;  cutY = 3.5;   nMeshPost = 150;
plotVelDecay = 0;   nDecay = 5;
plotPower = 0;

%filename
if ODE==0
    %filename = ['coneBubbleDeformable_x0=5_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUP) '_L=10_nStart=' num2str(nElemUP) '_rep=' num2str(repTypeUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    %filename = ['coneBubbleDeformable_ODE=' num2str(ODE) '_x0=5_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUP) '_L=10_nStart=' num2str(nElemUP) '_rep=' num2str(repTypeUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    filename = ['coneBubbleDeformable_massFlux=' num2str(massFlow) '_ODE=' num2str(ODE) '_x0=' num2str(x0upload) '_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUP) '_L=10_nStart=' num2str(nElemUP) '_rep=' num2str(repTypeUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
else
    filename = ['coneBubbleDeformable_ODE=' num2str(ODE) '_x0=' num2str(x0upload) '_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUP) '_L=10_nStart=' num2str(nElemUP) '_rep=' num2str(repTypeUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
end

%source
if res==1
   source = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==0
   source = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==2
   source = '~/Documents/MATLAB/test/coneWithDeformableBubble/resTest/';
end

BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorDeformablelBubble';
here = pwd;

saveDest = '/Users/Giacomo/Documents/research_notes/ETH_09052018/movies/frames/';
namePlot = 'plotConeExitFront';

%load file
load([source filename])
PARAM.ellipseShape = [0 0 0 0 1];
%PARAM.smoothingRep=1;

%number
if T(end)~=0
    loopUp = numel(T);
else
    loopUp = find(T==0,2);
    loopUp = loopUp(2)-1;
end

%loop
Ucone = zeros(numel(T),1);
manyPosMotor = zeros(numel(T),1);
Area = zeros(numel(T),1);
Volume = zeros(numel(T),1);
xcm = zeros(numel(T),1);
indT = find(T==0,2,'first');
if numel(indT)>1
    indT = indT(2);
    
    %clean data
    if ODE~=0
    count = 0;
    for i = 1:numel(T)
       if min(abs(T(i)-TsaveFirst))>1e-5
           i = i+1;
       else
           count = count+1;
       end
       Tnew(count) = T(i);
       Ynew{count} = Y{i};
    end
    T = Tnew';   Y = Ynew;
    loopUp = find(T==0,2);
    loopUp = loopUp(2)-1;
    end
    
    
else
    indT = numel(T);
end
Tplot = [Tplot nan];
tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp/180*pi linspace(0,1,100) linspace(-pi/2,pi/2,100)+thetaUp/180*pi linspace(0,pi,200)};
%compute power of cone with bubble with velocity 1
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
Udrag = 1;
PowerDragged = [];
Fdrag = [];
if plotPower==1
    [PowerDragged,Fdrag] = powerDraggedConeNoBubble(Udrag,tParametricPost(1:4),PARAM);
end
dXdt = zeros(loopUp,1);
PbubbleNoRep = zeros(loopUp,1);
PbubbleOnlyRep = zeros(loopUp,1);
PbubbleWithRep = zeros(loopUp,1);
Pmotor = zeros(loopUp,1);
countPlot = 1;
bubbleIsConfined = 0;
bubbleIsExiting = 0;
for i = 1:loopUp
    
   display([num2str(i) ' of ' num2str(indT)])
    
   %current shape and motor position
   Yhere = Y{i};
   xBubble = Yhere(1:2:end-2);
   yBubble = Yhere(2:2:end-1);
   posMotor = Yhere(end);
   
   %current location   TEMPORARY!!!
   PARAM.xStart(1:4) = PARAM.xStart(1:4)-PARAM.x0_Circle(2);
   PARAM.xEnd(1:4) = PARAM.xEnd(1:4)-PARAM.x0_Circle(2);
   PARAM.x0_Circle(1:4) = PARAM.x0_Circle(1:4)-PARAM.x0_Circle(2);
   PARAM_now = PARAM;
   PARAM_now.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
   PARAM_now.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
   PARAM_now.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;
   
   %compute curent shape
   [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);
   x{5} = xBubble';
   y{5} = yBubble';
   
   if bubbleIsConfined==0 && max(xBubble)<max(x{4})
           
           aBubble = max(xBubble)-min(xBubble);
           bBubble = 2*max(yBubble);
           Dbubble = (aBubble-bBubble)/(aBubble+bBubble);

           if Dbubble>stopValue*10
           
            indBubbleIsConfined = i+20;
            bubbleIsConfined = 1;
          
           end
           
   end
   
   if bubbleIsExiting==0 && max(xBubble)>max(x{4})
           
            indBubbleIsExiting = i;
            bubbleIsExiting = 1;
           
   end
   
   %motor velocity
   if ODE==0
    Vhere = V{i};
    Ucone(i) = Vhere(end);
    manyPosMotor(i) = Yhere(end);
   else
    manyPosMotor(i) = Yhere(end);
   end
   
   %bubble velocity
   uHere = Vhere(1:2:end-2); vHere = Vhere(2:2:end-1);
   xHere = Yhere(1:2:end-2); yHere = Yhere(2:2:end-1);
   Un = DropNormalVelocity(xHere',yHere',uHere,vHere);
   dXdt(i) = DropVelocityAxis(xHere',yHere',Un);
   
   %bubble surface area
   Area(i) = surf_gauss_vect(xBubble',yBubble');
   
   %bubble volume
   Volume(i) = axis_int_gauss_vect(xBubble',yBubble');
   
   %bubble center of mass
   xcm(i) = center_mass(xBubble',yBubble');
   
   if plotPower==1
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,CaUP,beta,PARAM);
       %PARAMvel = PARAM;
       %xVel = x;
       %yVel = y;
       PARAMvel.eliminateBlockPostprocessing = eliminateBlockPostprocessing;
       cd(here)
       
       PARAM.orderVariable = PARAM.orderVariableStokes;
       PARAM.orderGeometry = PARAM.orderGeometryStokes;
       PARAMvel.orderVariable = PARAMvel.orderVariableStokes;
       PARAMvel.orderGeometry = PARAMvel.orderGeometryStokes;
       
       %compute stress on cone from reoulsive forces on bubble
       [xCone,yCone] = getBlockCoordinates(xVel,yVel,PARAMvel,1);
       weightCone = integrationOnLineWeightAxis(xCone,yCone,PARAMvel.orderVariable(1),PARAMvel.orderGeometry(1),PARAMvel.SPlinesType(1));
       weightDrop = integrationOnLineWeightAxis(xHere',yHere',PARAMvel.orderVariable(5),PARAMvel.orderGeometry(5),PARAMvel.SPlinesType(5));
       
       %bubble
       [fBCx,fBCy,dfXrep,dfYrep] = computeStressOnBubble(xVel,yVel,5,PARAMvel);
       PbubbleNoRep(i) = weightDrop*(fBCx'.*uHere + fBCy'.*vHere);
       PbubbleOnlyRep(i) = weightDrop*(dfXrep.*uHere + dfYrep.*vHere);
       PbubbleWithRep(i) = weightDrop*((fBCx'+dfXrep).*uHere + (fBCy'+dfYrep).*vHere);
       PconeExternal(i) = weightDrop*((fBCx'+dfXrep).*Ucone(i));
       
       %cone
       fxCone = yStokes(1:2:2*numel(xCone)-3);
       %FonCone = -weightDrop*(dfXrep+fBCx');
       FonCone = -weightCone*fxCone;
       Pmotor(i) = -FonCone*Ucone(i);
       
   end
   
   %plot shape
   if (plotMotorUp==1 && sum(i==(1:plotIteUp:numel(T)))) || T(i)==Tplot(countPlot)
        
       if i==1 && savePlot==1
            width = 1201;
            height = 500;
            figure('pos',[50 100 width height])
       end
       
        %block coordinates
        if savePlot==1
            figure(2)
        else
            figure(1)
        end
        plotGeometryStokes(x,y,threeD,theta3D,color3D,transp,0,PARAM_now)
        axis([-5 15 -7 7])
        grid on
        hold off
        axis off
        %title(['t=' num2str(T(i))])
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
            print('-dpng','-loose','-r400',[saveDest namePlot sprintf('%03d',round((i-1)/plotIteUp)) '.png'])
        end
        
        if T(i)==Tplot(countPlot)
            
            figure(30)
            hold on
            plotGeometryStokesV2(x,y,threeD,theta3D,color3D,transp,0,-(countPlot-1)*5,PARAM_now)
            
            %axis([-5 15 -7 7])
            grid on
            hold off
            axis off
            %title(['t=' num2str(T(i))])
            drawnow
        
        end
        
   end
   
   if stopSphere==1 && max(xBubble)>max(x{4})
       
       aBubble = max(xBubble)-min(xBubble);
       bBubble = 2*max(yBubble);
       Dbubble = (aBubble-bBubble)/(aBubble+bBubble);
       
       if Dbubble<stopValue
          
           disp('Bubble is out and almost spherical')
           break
           
       end
       
   end
   
   %plot velocity field
   if plotVelField==1 && T(i)==Tplot(countPlot)
       
       %create grid
       Xgrid = linspace(-cutX,cutX,nMeshPost)+shift;
       Ygrid = linspace(0,cutY,nMeshPost);
       [Xgrid,Ygrid] = meshgrid(Xgrid,Ygrid);
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,CaUP,beta,PARAM);
       PARAMvel.eliminateBlockPostprocessing = eliminateBlockPostprocessing;
       cd(here)
       
       %compute velocity field
       [row,col] = size(Xgrid);
       uxGrid = zeros(row,col);
       uyGrid = zeros(row,col);
       pressure = zeros(row,col);
       for kkk = 1:row
            [Xgrid(kkk,:),Ygrid(kkk,:),uxGrid(kkk,:),uyGrid(kkk,:),pressure(kkk,:)] = computeVelPressFieldV2(Xgrid(kkk,:),Ygrid(kkk,:),xVel,yVel,yStokes,yStokes(end),0,PARAMvel,substitutePoint,coeffSub,noIN);
       end
    
       %[Xgrid,Ygrid,ux,uy] = computeVelPressField(Xgrid,Ygrid,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,substitutePoint,coeffSub,0);
       
       figure
       if streamslicePLot==1
            plotFieldAxisV2(x,y,Xgrid,Ygrid,uxGrid,uyGrid,1.5,0.2,[1 1],0,'parula','white',PARAMvel)
       else
           %countPlot = 2
           if countPlot==1
               
               startX = {5+0.6*cos(linspace(0,pi,4)) 6.2 []};
               startY = {0.6*sin(linspace(0,pi,4)) 3.5 []};
               
           elseif countPlot==2
               
               startX = {[10*ones(1,2) 10] -3*ones(1,3) []};
               startY = {[linspace(0,0.8,2) 3.5] linspace(0,2.9,3) []};
               
           elseif countPlot==3
               
               startX = {10+2*cos(linspace(0,pi/2,4)) -3*ones(1,3) 10.6};
               startY = {sin(linspace(0,pi/2,4)) linspace(0,2.9,3) 1.2};
               
           elseif countPlot==4
               
               startX = {10+2*cos(linspace(0,pi/2,3)) -3*ones(1,3) 10};
               startY = {2*sin(linspace(0,pi/2,3)) linspace(0,2.9,3) 3.5};
           
           end
           
           plotFieldAxisStreamlines(xVel,yVel,Xgrid,Ygrid,uxGrid,uyGrid,1.5,[1 1],0,'parula','white',startX,startY,yStokes,PARAMvel)
       
       end
       axis equal
       
       figure
       plotFieldAxisNoMap(x,y,Xgrid,Ygrid,uxGrid,uyGrid,1.5,0.2,[1 1],0,'black',PARAMvel)
       axis equal
       
   end
   
   %plot velocity field
   if plotVelFieldSlice==1 && T(i)==Tplot(countPlot)
       
       %create grid
       Ygrid = linspace(0,Ytop,nSlice);
       Xgrid = xSlice*ones(1,numel(Ygrid));
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,PARAM);
       PARAMvel.eliminateBlockPostprocessing = eliminateBlockPostprocessing;
       cd(here)
       
       %on the right of the motor
       [~,~,uxSlice,uySlice] = computeVelPressFieldV2(Xgrid,Ygrid,xVel,yVel,yStokes,yStokes(end),0,PARAMvel,0,coeffSub);
       
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
   if plotVelDecay==1 && T(i)==Tplot(countPlot)
       
%        %create grid
%        Xgrid = linspace(-cutX,cutX,nMeshPost)+shift;
%        Ygrid = linspace(0,cutY,nMeshPost);
%        [Xgrid,Ygrid] = meshgrid(Xgrid,Ygrid);
       
       %compute
       cd(BEM)
       [~,yStokes,xVel,yVel,PARAMvel] = computeVelocityDeformableBubbleODE(T(i),Y{i},tParametricBase,CaUP,beta,PARAM);
       PARAMvel.eliminateBlockPostprocessing = eliminateBlockPostprocessing;
       cd(here)
       
       %on the right of the motor
       Xright = logspace(1,nDecay,100);
       Yright = zeros(1,numel(Xright));
       [~,~,uxRight,uyRight] = computeVelPressFieldV2(Xright,Yright,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub,0);
       
       %on the left of the motor
       Xleft = -logspace(1,nDecay,100);
       Yleft = zeros(1,numel(Xleft));
       [~,~,uxLeft,uyLeft] = computeVelPressFieldV2(Xleft,Yleft,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub,0);
       
       %upward
       Yupward = logspace(1,nDecay,100);
       Xupward = zeros(1,numel(Yupward));
       [~,~,uxUpward,uyUpward] = computeVelPressFieldV2(Xupward,Yupward,xVel,yVel,yStokes,[yStokes(end) 0],0,PARAMvel,0,coeffSub,0);
       
       figure
       loglog(Xright,sqrt(uxRight.^2+uyRight.^2),'-x')
       hold on
       loglog(Yupward,sqrt(uxUpward.^2+uyUpward.^2),'-x')
       loglog(-Xleft,sqrt(uxLeft.^2+uyLeft.^2),'-x')
       grid on
       xyLabelTex('\rho','|U|')
       %xlabel('\rho')
       %ylabel('|U|')
       title('Velocity decay')
       drawnow
       legend('Right','Up','Left','Location','Best')
       
   end
   
   if T(i)==Tplot(countPlot)
       
       countPlot = countPlot+1;
       
   end
    
end

%derivative
D1 = finiteDifference1D(i,[2 0],1);

%area for a sphere
Rarea = nthroot(3/4/pi*Volume,3);
AreaSphere = 4*pi*Rarea.^2;
dVdt = (D1*Volume(1:i))./(D1*T(1:i));

%plot bubble volume and its variation
if plotBubbleVolume==1
    
    figure
    subplot(2,1,1)
    plot(T(1:i),Volume(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('V')
    title('Volume')
    drawnow
    
    subplot(2,1,2)
    %figure
    plot(T(1:i),dVdt(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('dV/dt')
    title('Volume variation')
    drawnow
      
end

%plot bubble surface area
if plotBubbleArea==1
    
    %compute dR/dt
    %dRdt = 1/4/pi/Rarea.^2.*dVdt;
    %dAdtSphere = 8*pi*Rarea.*dRdt;
    
    figure
    subplot(2,2,1)
    %figure
    plot(T(1:i),Area(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('A')
    title('Surface area')
    drawnow
    
    %figure
    subplot(2,2,2)
    plot(T(1:i),Area(1:i)-AreaSphere(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('\Delta A')
    title('Excess surface area')
    drawnow
    
    dAdt = (D1*Area(1:i))./(D1*T(1:i));
    dAdtSphere = (D1*AreaSphere(1:i))./(D1*T(1:i));
    
    subplot(2,2,3)
    %figure
    plot(T(1:i),dAdt(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('dA/dt')
    title('Surface area variation')
    drawnow
    
    subplot(2,2,4)
    %figure
    plot(T(1:i),(dAdt(1:i)-dAdtSphere(1:i)),'-k')
    hold on
    plot(T(1:i),dAdtSphere,'-r')
    grid on
    xlabel('t')
    ylabel('d (\Delta S)/dt')
    title('Excess surface area variation')
    legend('Due to deformation','Due to inflation','Location','Best')
    drawnow
      
end

%plot bubble center of mass
if plotBubbleXcm==1
    
    overlap1 = zeros(numel(Tplot)-1,1);
    overlap2 = zeros(numel(Tplot)-1,1);
    overlap3 = zeros(numel(Tplot)-1,1);
    overlap4 = zeros(numel(Tplot)-1,1);
    for l = 1:numel(Tplot)-1
        
        [~,ind] = min(abs(Tplot(l)-T));
        overlap1(l) = xcm(ind);
        overlap2(l) = manyPosMotor(ind);
        overlap3(l) = Ucone(ind);
        overlap4(l) = Area(ind)-AreaSphere(ind);
        
    end

    figure
    plot(T(1:i),xcm(1:i),'-k')
    hold on
    plot(Tplot(1:end-1),overlap1,'.r','MarkerSize',30)
    grid on
    xlabel('$t$','Interpreter','latex')
    title('Bubble position','Interpreter','latex')
    ylabel('$z_{b}$','Interpreter','latex')
    drawnow
    
    figure
    plot(T(1:i),manyPosMotor(1:i),'-k')
    hold on
    plot(Tplot(1:end-1),overlap2,'.r','MarkerSize',30)
    grid on
    xlabel('$t$','Interpreter','latex')
    title('Cone position','Interpreter','latex')
    ylabel('$z_{c}$','Interpreter','latex')
    drawnow
    
    figure(16)
    hold on
    %subplot(2,1,2)
    %figure
    %plot(T(1:i),dXdt(1:i),'-k')
    plot(T(1:i),dXdt(1:i))
    grid on
    xlabel('t')
    ylabel('U_b')
    title('Bubble velocity')
    drawnow
      
end

%plot motor velocity
if plotConeVel==1
    
    overlap3 = zeros(numel(Tplot)-1,1);
    for l = 1:numel(Tplot)-1
        
        [~,ind] = min(abs(Tplot(l)-T));
        overlap3(l) = Ucone(ind);
        
    end
    
    figure(17)
    hold on
    if ODE~=0
        Ucone = (D1*manyPosMotor(1:i))./(D1*T(1:i));
    end
       
    %figure
    plot(T(1:i),abs(Ucone(1:i)),'k')
    grid on
    hold on
    plot(Tplot(1:end-1),overlap3,'.r','MarkerSize',30)
    xyLabelTex('t','|\dot{z}_c|')
    %xlabel('t')
    %ylabel('U_M')
    title('Cone velocity')
    plot(T(indBubbleIsExiting),abs(Ucone(indBubbleIsExiting)),'.r','MarkerSize',30)
    plot(T(indBubbleIsConfined),abs(Ucone(indBubbleIsConfined)),'.g','MarkerSize',30)
    legend('Cone velocity','(t_2,|z_c^2|)','(t_1,|z_c^1|)','Location','Best')
    drawnow
    
    % fit power law
    figure
    loglog(T(indBubbleIsConfined:indBubbleIsExiting)-T(indBubbleIsConfined),abs(Ucone(indBubbleIsConfined:indBubbleIsExiting)-Ucone(indBubbleIsConfined)),'k')
    hold on
    fitPowerLaw = fit(T(indBubbleIsConfined+1:indBubbleIsExiting)-T(indBubbleIsConfined),abs(Ucone(indBubbleIsConfined+1:indBubbleIsExiting)-Ucone(indBubbleIsConfined)),'power1');
    xxx = linspace(5,T(indBubbleIsExiting)-T(indBubbleIsConfined),100);
    yyy = fitPowerLaw.a*xxx.^fitPowerLaw.b;
    loglog(xxx,yyy,'--k')
    grid on
    xyLabelTex('t-t_1','|\dot{z}_c-\dot{z}_c^1|')
    legend('From simulation',['Fitting ~t^{' num2str(fitPowerLaw.b) '}'],'Location','Best')
       
end

%plot power
if plotPower==1
    
    Pbubble = -dXdt*FonCone;
    
    if eliminateBlockPostprocessing(1)==0 && eliminateBlockPostprocessing(2)==0
    
        figure(18)
        plot(T(1:i),PbubbleNoRep(1:i))
        hold on
        grid on
        plot(T(1:i),PbubbleOnlyRep(1:i))
        plot(T(1:i),Pmotor(1:i))
        plot(T(1:i),dAdt(1:i),'--')
        plot(T(1:i),dAdt(1:i)-Pmotor(1:i)+Pbubble,'--')
        plot(T(1:i),PbubbleWithRep(1:i),'--')
        plot(T(1:i),PbubbleWithRep(1:i)+Pmotor(1:i),'--')
        plot(T(1:i),Pbubble,'--')
        %plot(T(1:i),PbubbleWithRep(1:i))
        %plot(T(1:i),PbubbleOnlyRep(1:i)+PbubbleNoRep(1:i),'--')
        legend('PbubbleNoRep','PbubbleOnlyRep','Pcone','dAdt','dAdt+Pcone+Pbubble','PbubbleWithRep','PdissTOT','Location','Best')
        xyLabelTex('t','P')
        title('Power')
    
    elseif eliminateBlockPostprocessing(1)==1 && eliminateBlockPostprocessing(2)==0
        
        %only bubble
        figure(18)
        plot(T(1:i),PbubbleNoRep(1:i))
        hold on
        grid on
        plot(T(1:i),dAdt(1:i),'x')
        plot(T(1:i),PbubbleOnlyRep(1:i))
        plot(T(1:i),Pbubble(1:i),'x')
        plot(T(1:i),PbubbleWithRep(1:i),'x')
        plot(T(1:i),PbubbleOnlyRep(1:i)+PbubbleNoRep(1:i),'o')
        legend('Pbubble','Pdiss','Location','Best')
        xyLabelTex('t','P')
        title('Power')
        
    elseif eliminateBlockPostprocessing(1)==0 && eliminateBlockPostprocessing(2)==1
        
        %only cone
        figure(18)
        plot(T(1:i),PconeExternal(1:i))
        hold on
        grid on
        plot(T(1:i),-Pmotor(1:i),'x')
        legend('Pcone','-Pdiss','Location','Best')
        xyLabelTex('t','P')
        title('Power')
    
    else
        
        error('Cannot eliminate all blocks')
        
    end

%     semilogy(T(1:i),PbubbleNoRep(1:i))
%     hold on
%     grid on
%     semilogy(T(1:i),PbubbleOnlyRep(1:i))
%     semilogy(T(1:i),Pmotor(1:i))
%     semilogy(T(1:i),dAdt(1:i),'--')
%     legend('PbubbleNoRep','PbubbleOnlyRep','Pcone','dAdt','Location','Best')

       
end

%check power dragged cone
DTunit = 1/Udrag;
EnergyPerUnitDisp = PowerDragged*DTunit;

%check power real cone
DT = T(i)-T(1);
DeltaZc = abs(manyPosMotor(i)-manyPosMotor(1));
deltaE = 99000/8.314/300;    %nondimensional energy jump for the reaction
Pchemical = deltaE*CaUP;
avgVel = DeltaZc/DT;
EnergySimulationPerUnitDisp = Pchemical*(1/avgVel);

efficiency = EnergyPerUnitDisp/EnergySimulationPerUnitDisp;
disp(['Efficiency eta=' num2str(efficiency)])

%check average velocity
Vavg = abs(manyPosMotor(i)/T(i));
disp(['Average cone velocity is Vavg=' num2str(Vavg)])

try
    disp(['Simulation time is ' num2str(simulationTime/60/60) ' hours'])
catch
    disp('Simulation is still running')
end













