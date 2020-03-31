%post processing spherical motor

clear variables
%close all

%parameters
loopUp = 2000;  
dtUp = 0.05;
res = 0;
thetaUp = pi/16;
minSizeElemRemeshUP = 1e-3/5;
coeffDistUP = 2;
remeshUP = 1;
remeshStepUP = 1;
rBubbleUP = 0.8;
xcmBubbleUP = 5;
repUP = 0;
nPerLenghtUP = 10;

%options
stopNow = loopUp-1;
plotMotorUp = 0;    plotIteUp = 10;
plotBubbleVel = 1;
plotMotorVel = 1;
plotElem = 1;
plotDistanceWallBubble = 0;
plotMesh = 0;

%shapshots
tSnapshot = 200;
%tSnapshot = 0:25:100;
plotVelField = 0;
plotVelDecay = 0;   distDecay = logspace(1,5,100);

%build grid
meshSize = 60;  substituitePoint = 0;
cutX = 10;  cutY = 10; shift=5;
x = linspace(-cutX,cutX,meshSize)+shift;
y = linspace(0,cutY,meshSize*round(cutY/cutX));
[Xgrid,Ygrid] = meshgrid(x,y);

if res==1
    source = '~/Documents/MATLAB/test/remeshGeometryBlocks/resTest/';
elseif res==0
    source = '~/Documents/MATLAB/droplet_simulations/server/';
end

%filename = ['testConeRemesh_theta=' num2str(thetaUp) '_loop=' num2str(loopUp) '_dt=' num2str(dtUp) '.mat'];
%filename = ['testConeRemesh_nPerL=' num2str(nPerLenghtUP) '_rep=' num2str(repUP) '_remesh=' num2str(remeshUP) '_step=' num2str(remeshStepUP) '_minSize=' num2str(minSizeElemRemeshUP) '_coeffDist=' num2str(coeffDistUP) '_rIN=' num2str(rBubbleUP) '_xcmIN=' num2str(xcmBubbleUP) '_theta=' num2str(thetaUp) '_loop=' num2str(loopUp) '_dt=' num2str(dtUp) '.mat'];
filename = ['testConeRemesh_rep=' num2str(repUP) '_remesh=' num2str(remeshUP) '_step=' num2str(remeshStepUP) '_minSize=' num2str(minSizeElemRemeshUP) '_coeffDist=' num2str(coeffDistUP) '_rIN=' num2str(rBubbleUP) '_xcmIN=' num2str(xcmBubbleUP) '_theta=' num2str(thetaUp) '_loop=' num2str(loopUp) '_dt=' num2str(dtUp) '.mat'];

%upload
load([source filename])

dist = zeros(stopNow,1);
distWallBubble = zeros(stopNow,1);
countSnapshots = 1;
for i = 1:loop
    
    display([num2str(i) ' of ' num2str(loop)])
    
    xNow = manyX{i};
    yNow = manyY{i};
    
    if i==stopNow
        display('STOP')
        break
    end
    
    if isempty(xNow)
        
        display('No data')
        break
        
    end
    
    if plotDistanceWallBubble==1

        %distance wall-bubble
        xCircle = xNow{5};
        PARAM.x0_Circle(5) = (max(xCircle)+min(xCircle))/2;
        PARAM.rArc = [0 0 0 0 (max(xCircle)-min(xCircle))/2];
        [dist1,~,distWallBubble1] = panelDistance(xNow,yNow,5,3,PARAM);
        [dist2,~,distWallBubble2] = panelDistance(xNow,yNow,5,3,PARAM);
        dist(i) = min(dist1,dist2);
        distWallBubble(i) = min(distWallBubble1,distWallBubble2);
    
    end
    
    if plotMotorUp==1 && sum(i==(1:plotIteUp:loop))
        
        %block coordinates
        figure(1)
        plotGeometryDrop(xNow,yNow,PARAM,plotMesh)
        axis([-2 12 -7 7])
        grid on
        drawnow
        hold off
        
        x3 = xNow{3};
        x4 = xNow{4};
        y3 = yNow{3};
        y4 = yNow{4};
        if abs(x3(end)-x4(1))>1e-5
            error('Bug in remesh')
        end
        
    end
    
    %plot snapshot
    if dt*(i-1)==tSnapshot(countSnapshots)
        
        if numel(tSnapshot)>1
        countSnapshots = countSnapshots+1;
        end
        
       figure(1)
       
       if countSnapshots>1
           hold on
       end
       
       if plotVelField==1
           
           PARAM.n(3:5) = [elemWall(i) elemArc(i) elemBubble(i)]-1;
           solution = BEM_Stokes(xNow,yNow,PARAM);
           [Xsing,Ysing,ux,uy,p] = computeVelPressField(Xgrid,Ygrid,xNow,yNow,solution,0,PARAM,yStokes(end-1),substituitePoint);
           streamslice(Xsing,Ysing,ux,uy)
           hold on
           contourf(Xsing,-Ysing,sqrt(ux.^2+uy.^2))
           xlabel('z')
           ylabel('r')
           colorbar
           title('Velocity field')
           
       end
       
       hold on
       plotGeometryDrop(xNow,yNow,PARAM,plotMesh)
       axis([-cutX+shift cutX+shift -cutY cutY])
       %grid on
       xlabel('z')
       ylabel('r')
       drawnow
       axis off
       
       if plotVelDecay==1
           
           if plotVelField==0
            PARAM.n(3:5) = [elemWall(i) elemArc(i) elemBubble(i)]-1;
            solution = BEM_Stokes(xNow,yNow,PARAM);
           end
           
           XdecayRight = distDecay;
           YdecayRight = zeros(numel(distDecay),1);
           [~,~,uxDecayFront,uyDecayFront,pDecayFront] = computeVelPressField(XdecayRight,YdecayRight,xNow,yNow,solution,0,PARAM,0,0);
           
           XdecayLeft = -distDecay;
           YdecayLeft = zeros(numel(distDecay),1);
           [~,~,uxDecayLeft,uyDecayLeft,pDecayLeft] = computeVelPressField(XdecayLeft,YdecayLeft,xNow,yNow,solution,0,PARAM,0,0);
           
           XdecaySide = zeros(numel(distDecay),1);
           YdecaySide = distDecay;
           [~,~,uxDecaySide,uyDecaySide,pDecaySide] = computeVelPressField(XdecaySide,YdecaySide,xNow,yNow,solution,0,PARAM,0,0);
           
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
    
end

%check if bubble has exited the pipe
xNow = manyX{i-1};
yNow = manyY{i-1};
xUP = xNow{1};
yUP = yNow{1};
xCircle = xNow{5};
yCircle = yNow{5};
xcm = (max(xCircle)+min(xCircle))/2;
[height,indMax] = max(yUP);
width = xUP(indMax);
C0 = height-tan(thetaCone)*width;
x0 = -C0/tan(thetaCone);
l = height/sin(thetaCone);
exitX = l/cos(thetaCone) + x0;
if xcm>exitX
    display('Bubble has exited the cone')
else
    display('Bubble has NOT exited the cone')
end
hold on
% plot(xcm,0,'or')
% plot(exitX,0,'xk')

if i==loop
        time = (0:i-1)*dt;
elseif i<loop
        time = (0:i-2)*dt;
end

if plotBubbleVel==1
    
    figure
    plot(time,manyUbubble(1:numel(time)))
    xlabel('t')
    ylabel('U_b')
    title('Bubble velocity')
    grid on
    
end

if plotMotorVel==1
    
    figure
    plot(time,manyUmotor(1:numel(time)))
    xlabel('t')
    ylabel('U_m')
    title('Motor velocity')
    grid on
    
end

if plotElem==1
    
    figure
    plot(time,elemWall(1:numel(time)))
    hold on
    plot(time,elemArc(1:numel(time)))
    plot(time,elemBubble(1:numel(time)))
    xlabel('t')
    ylabel('n')
    legend('straight wall','arc wall','bubble','Location','Best')
    title('Number of elements')
    grid on
    
end

if plotDistanceWallBubble==1
    
    figure
    %semilogy((0:i-2)*dt,distWallBubble)
    %hold on
    
    semilogy(time,dist(1:numel(time)),'-')
    xlabel('t')
    ylabel('d')
    %legend('min node distance','min shape dist','Location','Best')
    title('Distance wall-bubble')
    grid on
    
end











