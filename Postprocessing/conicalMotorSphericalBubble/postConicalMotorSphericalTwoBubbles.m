%post processing spherical motor

clear variables
close all

%parameters
TendUP = 1;
dtUp = 1e-3;
res = 1;
%thetaUp = linspace(pi/64,pi/16,100);
thetaUp = 2/180*pi;
coeffDistUP = 2;
remeshUP = 1;
repUP = 3;
coeffRepUP = 1e7;
distRepUP = 0.1;
nPerLenghtUP = 5;
tolUp = 1e-3;
ODEup = 2;
betaUp = 1;
Hup = 0.1;
Lup = 10;
BNup = 20;
nCycle = 10;
cycleType = 3;

%options
plotMotorUp = 1;    plotIteUp = 5;     savePlot = 1;
threeD = 1; theta3D = {[0 pi] [0 2*pi] [0 2*pi]};    color3D = {10 [1 1 1] [1 1 1]};   transp = [1 0.3 0.3];
plotBubbleDisp = 1;  plotBubbleVel = 0;
plotMotorDisp = 1;    plotMotorVel = 1;
plotElem = 0;
plotDistanceWallBubble = 1;
plotMesh = 0;
plotDOF = 0;
plotVolume = 0;
plotBubble1Radius = 0;

%shapshots
%tSnapshot = [0 .50 1 1.50 2]*0.2;
tSnapshot = 10;
plotVelField = 0;   plotPressureField = 0;  streamSLICE = 0;
plotPressureLine = 0;  plotStressBubble = 0;
plotVelDecay = 0;   distDecay = logspace(1,5,100);
plotConcField = 0;
plotConcLine = 0;

%build grid
meshSize = 45;  substitutePoint = 1;   coeffSub = 1;
cutX = 8;  cutY = 3; shift=6;
x = linspace(-cutX,cutX,meshSize)+shift;
y = linspace(0,cutY,meshSize);
[Xgrid,Ygrid] = meshgrid(x,y);

if res==1
    source = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==0
    source = '~/Documents/MATLAB/droplet_simulations/server/';
end

BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorWithDiffusionSphericalBubble';
here = pwd;

color = get(gca,'ColorOrder');
saveDest = '/Users/Giacomo/Documents/research_notes/micromotorWithDIffusion/movies/frames/';
%namePlot = 'plotConeVelocityOnly';

%destination for saving plot
%saveDest = '/Users/Giacomo/Documents/research_notes/bioLunch_Cambridge/movie/frame/';

if threeD==1
    namePlot = 'plotConeTwoBubbles';
else
    namePlot = 'plotCone';
end

if isempty(nCycle)
    
    filename = ['coneWithDiffusion_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUp) '_Hcc=' num2str(Hup) '_BN=' num2str(BNup) '_L=' num2str(Lup) '_nPerL=' num2str(nPerLenghtUP) '_rep=' num2str(repUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=' num2str(distRepUP) '_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];

    %upload
    load([source filename])

else
   
    filename = ['coneWithDiffusion_cycleType=' num2str(cycleType) '_nCycle=' num2str(nCycle) '_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUp) '_Hcc=' num2str(Hup) '_BN=' num2str(BNup) '_L=' num2str(Lup) '_nPerL=' num2str(nPerLenghtUP) '_rep=' num2str(repUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=' num2str(distRepUP) '_theta=' num2str(thetaUp) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    
    %upload
    load([source filename])
    
    %paste data
    for lll = 1:nCycle
        
        Tnow = manyT{lll};
        Ynow = manyY{lll};
        TinitCycle(lll) = Tnow(1);
        if lll==1
            Ynow = [Ynow(:,1:2) zeros(numel(Tnow),1) Ynow(:,3) zeros(numel(Tnow),1)];
            T = Tnow;
            Y = Ynow;
        else
            nBubble = (size(Ynow,2)-1)/2;
            
            if nBubble==2
                T = [T; Tnow(2:end)];
                Y = [Y; Ynow(2:end,:)];
            elseif nBubble==1
                Ynow = [Ynow(2:end,1:2) zeros(numel(Tnow)-1,1) Ynow(2:end,3) zeros(numel(Tnow)-1,1)];
                T = [T; Tnow(2:end)];
                Y = [Y; Ynow];
            else
                error('It is possible to have only one or two bubble')
            end
        end
        
    end
    
end
TinitCycle = [TinitCycle nan];

Vol1 = zeros(size(Y,1),1);
dist1 = zeros(size(Y,1),1);
distWallBubble1 = zeros(size(Y,1),1);
Vol2 = zeros(size(Y,1),1);
dist2 = zeros(size(Y,1),1);
distWallBubble2 = zeros(size(Y,1),1);
countSnapshots = 1;
bubbleBefore = 0;
countElim = 0;
countNucl = 0;
countSit = 1;
countDummy = 0;
for i = 1:numel(T)
    
    display([num2str(i) ' of ' num2str(numel(T))])
    
    %find our if there is one or two bubble at this instant
    if Y(i,3)==0 && Y(i,5)==0
        bubbleNow = 1;
        tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp linspace(0,1,100) linspace(-pi/2,pi/2,10)+thetaUp linspace(0,pi,100)};
        PARAM.n = [100 10 100 10 100];
        PARAM.panels = [4 1];
        
        if bubbleBefore==2
           
            countElim = countElim+1;
            Telim(countElim) = T(i);
            Yelim(countElim) = Y(i,1);
            countSit = countSit+1;
            
        end
        
    else
        bubbleNow = 2;
        tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp linspace(0,1,100) linspace(-pi/2,pi/2,10)+thetaUp linspace(0,pi,100) linspace(0,pi,100)};
        PARAM.n = [100 10 100 10 100 100];
        PARAM.panels = [4 1 1];
        
        if bubbleBefore==1
           
            countNucl = countNucl+1;
            Tnucl(countNucl) = T(i);
            Ynucl(countNucl) = Y(i,1);
            countSit = countSit+1;
            
        end
        
    end
    bubbleBefore = bubbleNow;
    
    %current location
    PARAM_now = PARAM;
    PARAM_now.xStart(1:4) = PARAM.xStart(1:4)+Y(i,1);
    PARAM_now.xEnd(1:4) = PARAM.xEnd(1:4)+Y(i,1);
    PARAM_now.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+Y(i,1);
    PARAM_now.x0_Circle(5) = Y(i,2);
    PARAM_now.rArc(5) = Y(i,4);
    if bubbleNow==2
        PARAM_now.x0_Circle(6) = Y(i,3);
        PARAM_now.rArc(6) = Y(i,5);
    end
    
    %compute curent shape
    [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);
    
    if plotDistanceWallBubble==1
        %compute rep forces
        PARAM_now2 = PARAM_now;
        PARAM_now2.orderVariable = PARAM_now2.orderVariableStokes;
        PARAM_now2.orderGeometry = PARAM_now2.orderGeometryStokes;
        [~,distWallBubble(i)] = computeRepulsiveForce2(x,y,1,numel(PARAM.panels),PARAM_now2);
        
    end
    
    if plotBubble1Radius==1 && cycleType==2
        
        if countSit==1
            bubbleRadius1(i) = Y(i,4);
            bubbleRadius2(i) = nan;
            bubbleRadius3(i) = nan;
        elseif countSit==2
            bubbleRadius1(i) = Y(i,4);
            bubbleRadius2(i) = Y(i,5);
            bubbleRadius3(i) = nan;
        elseif countSit==3
            bubbleRadius1(i) = nan;
            bubbleRadius2(i) = Y(i,4);
            bubbleRadius3(i) = nan;
        elseif countSit==4
            bubbleRadius1(i) = nan;
            bubbleRadius2(i) = Y(i,4);
            bubbleRadius3(i) = Y(i,5);
        elseif countSit==5
            bubbleRadius1(i) = nan;
            bubbleRadius2(i) = nan;
            bubbleRadius3(i) = Y(i,4);
        end
        
    elseif plotBubble1Radius==1 && cycleType==3
        
        
        if countDummy==1
            bubbleRadius1(i) = Y(i,4);
            bubbleRadius2(i) = nan;
            bubbleRadius3(i) = nan;
        elseif countDummy==2
            bubbleRadius1(i) = Y(i,4);
            bubbleRadius2(i) = Y(i,5);
            bubbleRadius3(i) = nan;
        elseif countDummy==3
            bubbleRadius1(i) = nan;
            bubbleRadius2(i) = Y(i,4);
            bubbleRadius3(i) = Y(i,5);
        elseif countDummy==0
            bubbleRadius1(i) = Y(i,4);
            bubbleRadius2(i) = nan;
            bubbleRadius3(i) = nan;
        end
        
        if TinitCycle(countDummy+1)==T(i)
            countDummy = countDummy+1;
        end
        
    end
    
    if plotMotorUp==1 && sum(i==(1:plotIteUp:size(Y,1)))
        
        %subplot(2,1,1)
        
        %block coordinates
        %figure(1)
        width = 1201;
        height = 500;
        if i==1
            figure('pos',[50 100 width height])
        end
        subplot(1,2,1)
        plotGeometry(x,y,threeD,theta3D,color3D,transp,PARAM_now)
        %axis([-2 12 -7 7])
        axis([-8 15 -5 5 -5 5])
        grid on
        hold off
        axis off
        %title(['\theta=' num2str(thetaUp)])
        drawnow
        
        %pause(0.1);
        
%         if i>1
%             %subplot(2,1,2)
%             figure(1)
%             D1 = finiteDifference1D(i,[2 0],1);
%             UmHere = D1*Y(1:i,1)./(D1*T(1:i));
%             plot(T(1:i),UmHere)
%             hold on
%             plot(T(i),UmHere(end),'.','MarkerSize',40,'Color',color(1,:))
%             grid on
%             xlabel('t')
%             ylabel('U_{cone}')
%             hold off
%             axis([0 2 -10 0])
%             drawnow
%         end

        if i>0
            subplot(1,2,2)
            %figure(1)
            plot(T(1:i),Y(1:i,1),'k')
            hold on
            plot(T(i),Y(i,1),'.k','MarkerSize',40)
            grid on
            xlabel('$t$','Interpreter','latex')
            ylabel('$z_c$','Interpreter','latex')
            hold off
            axis([0 2 -10 0])
            title('Cone displacement','Interpreter','latex')
            drawnow
        end
        
        if savePlot==1
            print('-dpng','-loose','-r400',[saveDest namePlot sprintf('%03d',round((i-1)/plotIteUp)) '.png'])
        end
        
    end
    
    %plot snapshot
    if abs(T(i)-tSnapshot(countSnapshots))<1e-5
        
        if countSnapshots>1
           hold on
        end
        
       figure(1)
       
       if plotVelField==1 || plotPressureField==1 || plotConcField==1
                      
           PARAM.SolveLaplace = 1;
           
           cd(BEM)
           if bubbleNow==1
                PARAMcomputeHere = PARAM;
                PARAMcomputeHere.remeshProximity = {[] [] 5 5 [3 4]};
                PARAMcomputeHere.blockType = [1 1];
                [~,yLaplace,yStokes,xVel,yVel,PARAMvel] = computeVelocityConicalMotorODE(T(i),Y(i,[1 2 4]),tParametricBase(1:5),nPerLenght,PARAMcomputeHere,Hcc,beta);
           elseif bubbleNow==2
                [~,yStokes,yLaplace,~,~,xVel,yVel,PARAMvel] = computeVelocityConicalMotorTwoBubblesODE(T(i),Y(i,:),tParametricBase,nPerLenght,PARAM,Hcc,beta);
           end
           cd(here)
           
           if threeD==0
               
               if plotVelField==1 || plotPressureField==1
               
               [Xsing,Ysing,ux,uy,p] = computeVelPressField(Xgrid,Ygrid,xVel,yVel,yStokes,yStokes(end-bubbleNow:end),0,PARAMvel,substitutePoint,coeffSub,1);
               
               if plotPressureField==0
                       [~,h1] = contourf(Xsing,Ysing,sqrt(ux.^2+uy.^2),500);
                       hold on
                       [~,h2] = contourf(Xsing,-Ysing,sqrt(ux.^2+uy.^2),500);
               elseif plotPressureField==1
                       [~,h1] = contourf(Xsing,Ysing,p,500);
                       hold on
                       [~,h2] = contourf(Xsing,-Ysing,p,500);
               end
                   
               if streamSLICE==1
                       hStream1 = streamslice(Xsing,Ysing,ux,uy,0.5);
                       hStream2 = streamslice(Xsing,-Ysing,ux,-uy,0.5);
               else
                       
                       %PUSHER
                       %startX1 = -100:33:100;  startY1 = 200*ones(1,7);
                       %startX2 = 20*ones(1,4);  startY2 = 0:3:10;
                       %startX3 = -10*ones(1,4);  startY3 = 0:3:10;
                       
                       %PULLER
                       %startX1 = 200*ones(1,6);  startY1 = 0:40:200;
                       %startX2 = -200*ones(1,6);  startY2 = 0:40:200;
                       %startX3 = 0:-0.1:-1;  startY3 = 0:2:20;
                       
                       %snap 1 T=0.155
%                        startX1 = -2*ones(1,3);  startY1 = 0:1:2.9;
%                        thetaHere = linspace(0,pi/4,4);
%                        startX2 = Y(i,2)+1.1*Y(i,4)*cos(thetaHere);  startY2 = 1.1*Y(i,4)*sin(thetaHere);
%                        startX3 = 4.2;  startY3 = 0;
                       
                       %snap 2 T=0.22
%                        thetaHere = linspace(0,pi,8);
%                        startX1 = Y(i,3)+1.1*Y(i,5)*cos(thetaHere);  startY1 = 1.1*Y(i,5)*sin(thetaHere);
%                        thetaHere = linspace(0,pi/2,4);
%                        startX2 = Y(i,2)+1.1*Y(i,4)*cos(thetaHere);  startY2 = 1.1*Y(i,4)*sin(thetaHere);
%                        startX3 = [];  startY3 = [];
                       
                       %snap 3 T=0.3
%                        thetaHere = linspace(0,pi/4,3);
%                        startX1 = Y(i,3)+1.2*Y(i,5)*cos(thetaHere);  startY1 = 1.2*Y(i,5)*sin(thetaHere);
%                        thetaHere = linspace(0,pi/2,4);
%                        startX2 = Y(i,2)+1.1*Y(i,4)*cos(thetaHere);  startY2 = 1.1*Y(i,4)*sin(thetaHere);
%                        startX3 = -2*ones(1,4);  startY3 = 0:1:3;
                       
                       %snap 4 T=0.31
                       startX1 = -2*ones(1,4);  startY1 = 0:1:3;
                       thetaHere = linspace(0,pi/4,4);
                       startX2 = Y(i,2)+1.1*Y(i,4)*cos(thetaHere);  startY2 = 1.1*Y(i,4)*sin(thetaHere);
                       startX3 = Y(i,3)+1.1*Y(i,5)*cos(thetaHere);  startY3 = 1.1*Y(i,5)*sin(thetaHere);
                       %startX3 = [];  startY3 = [];
                       
                       %Close look, unconfined bubble
%                        startX1 = [];  startY1 = [];
%                        thetaHere = linspace(0,pi,6);
%                        rForPlot = Y(i,3)+0.1;
%                        startX2 = Y(i,2)+rForPlot*cos(thetaHere);  startY2 = rForPlot*sin(thetaHere);
%                        startX3 = [];  startY3 = [];
                       
                       optionStream = [0.1 1000];
                       %optionStream = [];
                       hStream1 = streamline(Xsing,Ysing,ux,uy,[startX1 startX2 startX3],[startY1 startY2 startY3],optionStream);
                       hStream2 = streamline(Xsing,-Ysing,ux,-uy,[startX1 startX2 startX3],-[startY1 startY2 startY3]);
                       xyCoord = stream2(Xsing,Ysing,ux,uy,[startX1 startX2 startX3],[startY1 startY2 startY3]);
                       %plot(startX3,startY3,'or')
               end
                   %hStream2 = streamline(Xsing,-Ysing,ux,-uy,[200*ones(1,6) -200*ones(1,6)],[0:-40:-200 0:-40:-200]);
                   xlabel('z')
                   ylabel('r')
                   colorbar
               if plotPressureField==1
                       title(['Pressure field t=' num2str(tSnapshot(countSnapshots))])
               else
                       title(['Velocity field t=' num2str(tSnapshot(countSnapshots))])
               end

                   set(hStream1,'LineWidth',1.5,'Color','white')
                   set(hStream2,'LineWidth',1.5,'Color','white')
                   set(h1,'LineColor','none')
                   set(h2,'LineColor','none')
                   caxis([0 400])
                   hold on
                   
               if streamSLICE==0
                   axis(3/2*[-2 12 -3 3])
                   %add arrows to stremlines
                   for zzz = 1:numel(xyCoord)
                       
                      xyCoordHere = xyCoord{zzz};
                      %xCoord = xyCoordHere(50:240:end,1);
                      %yCoord = xyCoordHere(50:240:end,2);
                      sizeHere = size(xyCoordHere);
                      xCoord = xyCoordHere(ceil(sizeHere(1)/3),1);
                      yCoord = xyCoordHere(ceil(sizeHere(1)/3),2);
                      
                      %compute velocity
                      [~,~,uxLine,uyLine,pLine] = computeVelPressField(xCoord,yCoord,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub,1);
                      
                      %quiver(xCoord,yCoord,uxLine,uyLine,'r','LineWidth',3)
                      %scaleHere = 5;
                      %quiver(xCoord,yCoord,scaleHere*uxLine./sqrt(uxLine.^2+uyLine.^2),scaleHere*uyLine./sqrt(uxLine.^2+uyLine.^2),'r','LineWidth',3,'AutoScale','Off','MaxHeadSize',10)
                      
                      dx = uxLine./sqrt(uxLine.^2+uyLine.^2)/100;
                      dy = uyLine./sqrt(uxLine.^2+uyLine.^2)/100;
                      dlHere = sqrt(diff(xCoord).^2+diff(yCoord).^2);
                      distOrigin = sqrt(xCoord.^2+yCoord.^2);
                      for sss = 1:numel(xCoord)
                          
                          if sss==1 && distOrigin(sss)>0
                              arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','w')
                              arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','w')
                          end
                          
                          if sss>1
                              if dlHere(sss-1)>40 && distOrigin(sss)>0
                                arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','r')
                                arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','r')
                              end
                          end
                      
                      end
                       
                   end
               end
               axis([min(min(Xgrid)) max(max(Xgrid)) min(min(Ygrid)) max(max(Ygrid))])
               
               elseif plotConcField==1
                  
                   %[Xsing,Ysing,ux,uy,p] = computeConcentrationField(Xgrid,Ygrid,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,substitutePoint,coeffSub);
                   [Xsing,Ysing,PHIfield] = computeConcentrationField(Xgrid,Ygrid,xVel,yVel,yLaplace,PARAMvel,substitutePoint,coeffSub);
                   
                   xlabel('z')
                   ylabel('r')
                   colorbar
                   title(['Concentration field t=' num2str(tSnapshot(countSnapshots))])
                   [~,h1] = contourf(Xsing,Ysing,PHIfield,500);
                   hold on
                   [~,h2] = contourf(Xsing,-Ysing,PHIfield,500);
                   set(h1,'LineColor','none')
                   set(h2,'LineColor','none')
                   colorbar
                   caxis([0 20])
                   hold on
                   
               end
           
           else
               
               %build 3D grid
               zGrid = linspace(-cutX,cutX,meshSize)+shift;
               yGrid = linspace(-cutY,cutY,round(meshSize*cutY/cutX));
               xGrid = linspace(-cutY,cutY,round(meshSize*cutY/cutX)+1);
               [Zgrid,Xgrid,Ygrid] = meshgrid(zGrid,xGrid,yGrid);
               
               %allocation
               uz = zeros(numel(xGrid),numel(zGrid),numel(yGrid));
               ur = zeros(numel(xGrid),numel(zGrid),numel(yGrid));
               theta = zeros(numel(xGrid),numel(yGrid),numel(yGrid));
               
               %compute on 3D grid
               for iii = 1:numel(xGrid)
                   for kkk = 1:numel(yGrid)
                       for lll = 1:numel(zGrid)
                           zHere = Zgrid(iii,lll,kkk);
                           rHere = sqrt(Xgrid(iii,lll,kkk)^2+Ygrid(iii,lll,kkk));
                           theta(iii,lll,kkk) = atan(Ygrid(iii,lll,kkk)./Xgrid(iii,lll,kkk)) + pi*((Xgrid(iii,lll,kkk)<0)&&(Ygrid(iii,lll,kkk)<0)) + pi*((Xgrid(iii,lll,kkk)<0)&&(Ygrid(iii,lll,kkk)>0));
                            [~,~,uz(iii,lll,kkk),ur(iii,lll,kkk),p] = computeVelPressField(zHere,rHere,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
                            if isnan(theta(iii,lll,kkk))
                                theta = 0;
                            end
                       end
                   end
               end
               
               %transform data into 3D
               ux = ur.*cos(theta);
               uy = ur.*sin(theta);
               
               %sx = linspace(min(xGrid),max(xGrid),5);
               %sy = linspace(min(xGrid),max(xGrid),5);
               thetaPlot = linspace(0,2*pi,6);
               sx = 2*cos(thetaPlot);
               sy = 2*sin(thetaPlot);
               sz = -4;
               [sz,sx,sy] = meshgrid(sz,sx,sy);
               %streamtube(Zgrid,Xgrid,Ygrid,uz,ux,uy,sz,sx,sy)
               streamtube(Zgrid,Xgrid,Ygrid,uz-yStokes(end-1),ux,uy,sz,sx,sy,[0.3 20])
               hold on
               shading interp
               camlight; lighting gouraud
               title(['Velocity field cone frame t=' num2str(tSnapshot(countSnapshots))])
               
           end
           
       end
       
       %hold on
       plotGeometry(x,y,threeD,theta3D,color3D,transp,PARAM_now);
       %axis([-cutX+shift cutX+shift -cutY cutY])
       %grid on
       xlabel('z')
       ylabel('r')
       %axis off
       if plotVelField==0
           title(['t=' num2str(tSnapshot(countSnapshots))])
       end
       drawnow
       
       if numel(tSnapshot)>1
            countSnapshots = countSnapshots+1;
       end
       
       if plotPressureLine==1 || plotStressBubble==1 || plotConcLine==1
           
           PARAM.SolveLaplace = 1;
           
           cd(BEM)
                [~,yLaplace,yStokes,xVel,yVel,PARAMvel] = computeVelocityConicalMotorODE(T(i),Y(i,:),tParametricBase,nPerLenght,PARAM,Hcc,beta);
           cd(here)
           
           %compute pressure
           if plotPressureLine==1
               
                Xline = linspace(-5,15,200);
                Yline = zeros(1,numel(Xline))+1e-2;
                [~,~,~,~,p] = computeVelPressField(Xline,Yline,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
            
                figure
                plot(Xline,p)
                grid on
                xlabel('z')
                ylabel('p')
                title('Pressure along the axis')
                
                thetaHere = linspace(0,pi,200);
                %xBubble = x{5}; yBubble = y{5};
                
                rHere = PARAMvel.rArc(5) + 0.05;
                %maxElemHere = (pi*rHere)/PARAMvel.n(5);
                %if maxElemHere>
                %rHere = rHere + ;
                Xline = rHere*cos(thetaHere) + PARAMvel.x0_Circle(5);
                Yline = rHere*sin(thetaHere);
                [~,~,~,~,p] = computeVelPressField(Xline,Yline,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
            
                figure(3)
                
                hold on
                plot(Xline,p)
                grid on
                xlabel('z')
                ylabel('p')
                title('Pressure close to the bubble')
                
                %compute integral
                pAvg = 0.5*trapz(thetaHere,p.*sin(thetaHere));
                
                %compute new beta
                %pAvg = 0;
                newBeta = (1e5+pAvg*1e-1)*1e-6/7.2/1e-2;
            
           end
           
           %compute stress on bubble
           if plotStressBubble==1
               
               figure(4)
               hold on
               fN = yStokes(2*sum(PARAMvel.n(1:4))+1:2:end-3);
               xVel = xVel{5};
               yVel = yVel{5};
               plot((xVel(1:end-1)+xVel(2:end))/2,fN)
               grid on
               xlabel('z')
               ylabel('f_n')
               title('Normal stress on bubble')
               drawnow
               
               %compute integral
               dl = sqrt(diff(xVel).^2+diff(yVel).^2);
               rInt = (yVel(1:end-1)+yVel(2:end))/2;
               intFN = 2*pi*sum(rInt.*dl.*(fN)');
               
           end
           
           if plotConcLine==1
               
               Xline = linspace(-5,15,100);
               Yline = zeros(1,numel(Xline));
               [~,~,PHIfield] = computeConcentrationField(Xline,Yline,xVel,yVel,yLaplace,PARAMvel,0,1);
               
               figure(5)
               hold on
               plot(Xline,PHIfield)
               grid on
               xlabel('z')
               ylabel('c')
               title('Concentration along the axis')
               drawnow
               
           end
           
       end
       
       if plotVelDecay==1
           
           if plotVelField==0
            [~,yLaplace,yStokes,xVel,yVel,PARAMvel] = computeVelocityConicalMotorODE(T(i),Y(i,:),tParametricBase,nPerLenght,PARAM,Hcc,beta);
           end
           
           XdecayRight = distDecay;
           YdecayRight = zeros(numel(distDecay),1);
           [~,~,uxDecayFront,uyDecayFront,pDecayFront] = computeVelPressField(XdecayRight,YdecayRight,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
           
           XdecayLeft = -distDecay;
           YdecayLeft = zeros(numel(distDecay),1);
           [~,~,uxDecayLeft,uyDecayLeft,pDecayLeft] = computeVelPressField(XdecayLeft,YdecayLeft,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
           
           XdecaySide = zeros(numel(distDecay),1);
           YdecaySide = distDecay;
           [~,~,uxDecaySide,uyDecaySide,pDecaySide] = computeVelPressField(XdecaySide,YdecaySide,xVel,yVel,yStokes,[yStokes(end-1) yStokes(end)],0,PARAMvel,0,coeffSub);
           
           figure
           loglog(XdecayRight,sqrt(uxDecayFront.^2+uyDecayFront.^2),'o-')
           hold on
           loglog(YdecaySide,sqrt(uxDecaySide.^2+uyDecaySide.^2),'o-')
           loglog(abs(XdecayLeft),sqrt(uxDecayLeft.^2+uyDecayLeft.^2),'o-')
           grid on
           xlabel('\rho')
           ylabel('|U|')
           legend('Right','Radial','Left','Location','Best')
           
       end
        
    end
    
    %compute bubble volume
    Vol(i) = 4/3*pi*PARAM_now.rArc(5)^3;
    
end

D1 = finiteDifference1D(numel(T),[2 0],1);

if plotBubbleDisp==1
    
    figure
    plot(T(1:i),Y(:,2))
    xlabel('t')
    ylabel('z_b')
    title('Bubble displacement')
    grid on
    
end

if plotMotorDisp==1
    
    figure
    plot(T(1:i),Y(:,1))
    xlabel('t')
    ylabel('z_m')
    title('Motor displacement')
    grid on
    
%     hold on
%     plot(Tnucl,Ynucl,'x')
%     plot(Telim,Yelim,'o')
%     ind = find(T==tSnapshot,2,'first');
%     %plot(tSnapshot,Y(ind,1),'s')
%     %legend('Disp','Nucl','Elim','Snapshot')
%     legend('Disp','Nucl','Elim')
    
end

if plotBubbleVel==1
    
    Ububble = D1*Y(:,2)./(D1*T);
    
    figure
    plot(T(1:i),Ububble(1:i))
    xlabel('t')
    ylabel('U_b')
    title('Bubble velocity')
    grid on
    
end

if plotMotorVel==1
    
    Umotor = D1*Y(:,1)./(D1*T);
    
    figure
    plot(T(1:i),Umotor(1:i))
    xlabel('t')
    ylabel('U_m')
    title('Motor velocity')
    grid on
    
end

if plotDistanceWallBubble==1
    
    figure
    %semilogy((0:i-2)*dt,distWallBubble)
    %hold on
    
    semilogy(T(1:i),distWallBubble(1:i),'-')
    xlabel('t')
    ylabel('d')
    %legend('min node distance','min shape dist','Location','Best')
    title('Distance wall-bubble')
    grid on
    
end

if plotDOF==1
    
   figure
   plot(T(1:i),Y(1:i))
   grid on
   xlabel('t')
   title('Degrees of freedom')
    
end

%plot volume
if plotVolume==1
        
       figure
       semilogy(T(1:i),Vol(1:i))
       xlabel('t')
       ylabel('Volume')
       title('Bubble volume')
       grid on
        
end

if plotBubble1Radius==1
        
    figure
    plot(T(1:i),bubbleRadius1(1:i));
    hold on
    plot(T(1:i),bubbleRadius2(1:i));
    plot(T(1:i),bubbleRadius3(1:i));
    grid on
    xlabel('t')
    ylabel('r_b')
        
end

%vAvg = 1/(T(end)-T(1))*trapz(T,Umotor);
vAvg = Y(i,1)/(T(i)-T(1));
disp(['Average motor velocity is ' num2str(abs(vAvg))])
disp(['Simulation time is ' num2str(simulationTime/60) ' min'])











