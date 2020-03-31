%plot DNS at some desired times

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

close all
clear variables

%problem parameters
%source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/';
source = '~/Documents/MATLAB/droplet_simulations/results/';
%source = '~/Documents/MATLAB/droplet_simulations/server/';
dest = '~/Documents/research_notes/micromotor/slideMeeting250417/movies/frames/';
theta = -0.02;
Ca = 0.01;
dt = 0.05;
visc = 0;
inflate = 0;
element = 82;
MotorFree = 0;
TypeRep = 6;
IntensRep = -1e2;
RepLenght = 2e-1;
Resize = 0;
alpha = 1.2;
L = 10;
Inpos = 5;
deflationWall = 1;

%upload data
%if isempty(deflationWall)
    %filename = [source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(TypeRep) '_IntenRep=' num2str(IntensRep) '_RepLenght=' num2str(RepLenght) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(element) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=' num2str(L) '_alpha=' num2str(alpha) '_InPos=' num2str(Inpos) '_RK=2.mat'];
%else
    filename = [source 'ConicalMotor_deflationWall=' num2str(deflationWall) '_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(TypeRep) '_IntenRep=' num2str(IntensRep) '_RepLenght=' num2str(RepLenght) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(element) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=' num2str(L) '_alpha=' num2str(alpha) '_InPos=' num2str(Inpos) '_RK=2.mat'];
%end
load(filename)
PARAM.BEM = '~/Documents/MATLAB/droplet_simulations/conicalMotor/';

%if isempty(deflationWall)
 %   PARAM.deflationWall = 0;
%end

%plotting options
step = 1;                                                      % plot every few iteration
PlotSnapshot = 1;   plotInitialMotor = 0;                       % plot snapshots
SaveSnapshot = 0;                                               % save snapshots
SubPlot = 0;                                                    % plot with subplot
SubPlotStress = 0;                                              % plot with subplot for stresses
PlotVelField = 0;   decompose = 0;                              % plot velocity field
SaveVelField = 0;                                               % save vel field
computePressureField = 0;
ComputeDropVel = 1;                                             % compute velocity of the droplet from the solution
ComputeDropPos = 0;                                             % compute droplet position
ComputeDropSurface = 1;                                         % compute surface area of the bubble
ComputeSurfaceVariation = 1;
ComputeDropVolume = 0;                                          % compute volume of the bubble
PlotDropVel = 0;                                                % plot velocity of the droplet from the solution
PlotStresses = 0;                                               % plot stresses distribution on the motor wall
PlotVelNow = 0;                                                 % plot motor velocity at this instant
SavePlotVelNow = 0;                                             % save motor velocity at this instant
PlotDisplaceNow = 0;                                            % plot motor velocity at this instant
SavePlotDisplaceNow = 0;                                        % save motor velocity at this instant
PlotVelSlice = 0;   xSlice = 8;     scale = 0; decreaseSliceY = 0.2;% plot vel fiedl only in one location
PlotPressureSlice = 0;   ySlice = 0; xLimit = [0 10];   thetaSlice = 0;           % plot vel fiedl only in one location
Zoom = 0;   xLim = [0.45 10];  yLim = [0.5 1.4];  MeshZ = 15;   % plot vel field zomming in one location
CheckDecay = 1;     nDecay = 100;   dist = 5;                   % check velocity decay
ModifyStresses = 0;     signPressure = -1;                      % modify stresses by adding arbitrary pressure
AddArtificialStokeslet = 0;                                     % balance net contribution form drop
checkCapillaryTwoThird = 0;     scaling = 2/3;  mult = 1e-1;    % check bretherton scaling
plotFinalCurvature = 0;                                         % plot curvature of the final shape
plotStresses = 1;

%options physics
FrameLab = 1;                                           % 0 plot in motor frame, 1 lab frame or 2 bubble frame
MASK = 0;                                               % eventually don't plot velocity inside the drop
MASKzoom = 0;                                           % eventually don't plot velocity inside the drop

%mesh field option
cutX = 10;   cutY = 10;   shift = 5;
xRange1 = -cutX + shift;   xRange2 = cutX + shift;
yRange1 = 0.0;   yRange2 = cutY;
MeshFine = 100/cutX;

%mesh field option
MeshFineSlice = 200;

%which time to plot
time = 590;
%time = linspace(1,1900,4);
%time = [1 3000 5250 7000];

%figure in pixel
width = 1400;
height = 600;

%first coordinate of motor
m = find(risa(2:end,1)==risa(1,1));
indNode = find(risb(:,1)==0,2,'first');
aMotorFirst = risa(1:m+1,1);
bMotorFirst = risb(1:m+1,1);

PARAM.computePressure = computePressureField;

%velocity data
ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1)-1;
v = zeros(1,ite/step);
displace = zeros(1,ite/step);
vDrop = zeros(1,ite/step);
dropPos = zeros(1,ite/step);
AreaDrop = zeros(1,ite/step);
xcmDrop = zeros(1,ite/step);
hBre = zeros(1,ite/step);
volumeDrop = zeros(1,ite/step);

for k = 1:ite/step
    %display(i)
    i = k*step+1;
    indVel = find(risy(:,i)==0,2,'first');
    v(k) = risy(indVel(2)-1,i);
    displace(k) = risa(1,i)-risa(1,1);
    
    if ComputeDropVel==1||ComputeDropPos==1
        
        display([num2str(i) ' of ' num2str(ite)])
        
        %often used
        m = find(risa(2:end,i)==risa(1,i));
        indNode = find(risb(:,i)==0,2,'first');
        q = indNode(2)-m-3;
        aWall = risa(1:m+1,i);          bWall = risb(1:m+1,i);
        aDrop = risa(m+2:m+q+2,i);      bDrop = risb(m+2:m+q+2,i);
        
        if ComputeDropVel==1
            VelInterfaceX = risy(2*m+1:2:2*(m+q)+1,i);
            VelInterfaceY = risy(2*m+2:2:2*(m+q)+2,i);

            %compute vector normal to the interface
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(aDrop',bDrop');
            N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
                -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

            %velocity normal to interface
            Unormal = VelInterfaceX'.*N(1,:) + VelInterfaceY'.*N(2,:);
            
            %compute velocity
            vDrop(k) = DropVelocityAxis(aDrop',bDrop',Unormal');
        
        end
        
        if ComputeDropPos==1
            
            dropPos(k) = center_mass(aDrop,bDrop);
            
        end
        
        if ComputeDropSurface==1
            
            %compute surface area
            AreaDrop(k) = surf_gauss_vect(aDrop',bDrop');
            %xcmDrop(k) = center_mass_gauss(aDrop',bDrop');
        
        end
        
        if ComputeDropVolume == 1
            
            %compute volume
            volumeDrop(k) = axis_int_gauss_vect(aDrop',bDrop');
            
        end
        
        if checkCapillaryTwoThird==1
            
            hBre(k) = min(sqrt((aDrop(round(q/2))-aWall).^2+(bDrop(round(q/2))-bWall).^2));
            
        end
        
%         figure(10)
%         plot(aDrop,bDrop)
%         axis equal
%         grid on
%         axis([0 5 -2 2])
%         drawnow
        
    end
    
end
tttHere = 0:PARAM.checkpoint*PARAM.deltaT*step:time(end)-PARAM.checkpoint*PARAM.deltaT;

if SubPlot==1
    figure(1)
    subplot((numel(time)+2)/3,3,1)
else
    figure(2)
end
plot(tttHere,v,'LineWidth',2)
%plot average velocity
hold on
vavg = sum(v)/numel(v);
plot(tttHere,vavg*ones(numel(v),1),'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
%axis([0 time(end) min(risy(m+q+2,:)) max(risy(m+q+2,:))])
ylabel('v')

if PlotSnapshot==1
    for i = 1:numel(time)

        display([num2str(i) ' of ' num2str(numel(time))])

        if i==1
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1);
        else
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1)-1;
        end

        %often used
        %m = PARAM.m;    %q = nbr_el(ite);
        m = find(risa(2:end,ite)==risa(1,ite));
        indNode = find(risb(:,ite)==0,2,'first');
        q = indNode(2)-m-3;

        aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
        aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);

        %velocity data
        if numel(time)==1
            vnow = v(1:ite/step-1);
            displaceNow = displace(1:ite/step-1);
        else
            vnow = v(1:ceil(ite/step));
            displaceNow = displace(1:ceil(ite/step));
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(i);
        
        if PlotVelNow==1
            
            figure(3)
            if numel(ttt)==numel(vnow)
                plot(ttt,vnow)
            else
                plot(ttt(2:end),vnow)
            end
            hold on
            plot(ttt(end),vnow(end),'.k','MarkerSize',40)
            hold off
            grid on
            xlabel('t')
            ylabel('v')
            title('Motor velocity')
            axis([0 time(end) -1e-3 1e-3])
            drawnow
            
            if SavePlotVelNow==1
                %fig = figure(3);
                %fig = gcf;
                %fig.Position = [400 200 width height];
                name = 'MotorVelocity';
                print('-dpng','-loose','-r100',[dest name sprintf('%04d',i) '.png'])
            end
            
        end
        
        if PlotDisplaceNow==1
            
            figure(4)
            if numel(ttt)==numel(displaceNow)
                plot(ttt,displaceNow)
            else
                plot(ttt(2:end),displaceNow)
            end
            hold on
            plot(ttt(end),displaceNow(end),'.k','MarkerSize',40)
            hold off
            grid on
            xlabel('t')
            ylabel('v')
            title('Motor displacement')
            axis([0 time(end) 0 0.35])
            drawnow
            
            if SavePlotDisplaceNow==1
                %fig = figure(3);
                %fig = gcf;
                %fig.Position = [400 200 width height];
                name = 'MotorDisplace';
                print('-dpng','-loose','-r100',[dest name sprintf('%04d',i) '.png'])
            end
            
        end
        
        figure(1)
        fig = figure(1);
        
        if SubPlot==1
            fig.Position = [400 200 width height];
            subplot((numel(time)+2)/3,3,i+2)
            %subplot((numel(time)+1)/2,2,i+1)
        else
            %figure(1)
            figure(1)
            hold on
            %hold off
        end
        
        if PlotVelField==1
            
            %grid for computing velocity
            x = linspace(xRange1,xRange2,MeshFine*(xRange2-xRange1));
            y = linspace(yRange1,yRange2,MeshFine*yRange2);
            [X,Y] = meshgrid(x,y);
            
            %select solution
            numSing = numel(aDrop)+numel(aMotor)-1;
            solution = risy(1:2*numSing+1,ite);
            
            %compute velocity field
            [U,V,uFrame,p] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',visc,X,Y,PARAM,FrameLab,MASK,decompose);
            
            streamslice(X,-Y,(U-uFrame),-V,'b')
            hold on
            contourf(X,-Y,sqrt((U-uFrame).^2+V.^2))
            colorbar
            streamslice(X,Y,U-uFrame,V,'b')
            if FrameLab==1
                title('Lab Frame')
            elseif FrameLab==0
                title('Motor Frame')
            elseif FrameLab==2
                title('Drop frame')
            end
            
        elseif PlotVelSlice==1
           
            %choose grid
            %xSlice = 10;
            %yRangeSlice1 = max(bMotor)-1;
            yRangeSlice1 = 0;
            yRangeSlice2 = max(bMotor)-decreaseSliceY;
            %yRangeSlice2 = 50;
            x = xSlice;
            y = linspace(yRangeSlice1,yRangeSlice2,MeshFineSlice*yRangeSlice2);
            [X,Y] = meshgrid(x,y);
            
            %select solution
            numSing = numel(aDrop)+numel(aMotor)-1;
            solution = risy(1:2*numSing+1,ite);
            
            %compute velocity field
            [U,V,uFrame] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',visc,X,Y,PARAM,FrameLab,MASK,decompose);
            
            %compute flow rate
            Qslice = 2*pi*trapz(Y,Y.*U);
            display(Qslice)
            
            if ComputeDropVel==1
                %flow rate assocaited to the bubble
                Qdrop = pi*vDrop(end);
                display(Qdrop)
            end
            
            if scale==1
                U = 2*U/max(abs(U));
                uFrame = 2*uFrame/max(abs(U));
            end
            
            plot(U+xSlice-uFrame,y,'b-')
            %quiver(X,-Y,(U-uFrame),-V,'b')
            hold on
            plot(U+xSlice-uFrame,-y,'b-')
            plot(X,Y,'m-')
            plot(X,-Y,'m-')
            %streamslice(X,Y,U-uFrame,V,'b')
            
        end
           
        plot(aDrop,bDrop,'r',aDrop,-bDrop,'r')
        xcm = center_mass(aDrop,bDrop);
        hold on
        if plotInitialMotor==1
            grey = [0.4 , 0.4 , 0.4] ;
            p = plot(aMotorFirst,bMotorFirst,aMotorFirst,-bMotorFirst);
            set (p , 'Color', grey)
        end
        plot(aMotor,bMotor,'k',aMotor,-bMotor,'k')
        %plot(xcm,0,'k.','MarkerSize',20)
        grid on
        axis equal
        if PlotVelField==1
            axis([xRange1 xRange2 -yRange2 yRange2])
        elseif PlotVelSlice==1
            %axis([xSlice-1 xSlice+1 -yRangeSlice2-0.5 yRangeSlice2+0.5])
            axis([-2 12 -4 4])
        else
            axis([-3 L+2 -3 3])
        end
        xlabel('x')
        ylabel('r')
        drawnow
        hold off
        
        
        if PlotVelSlice==1
            figure(3)
            plot(Y,U-uFrame)
            grid on
            xlabel('r')
            ylabel('V')
            if scale==1
                ylabel('V/V_m')
            end
        end
        
        if SaveSnapshot==1

            fig = gcf;
            axis off
            %fig.Position = [400 200 width height];
            name = 'MotorSnapshot';
            print('-dpng','-loose','-r100',[dest name sprintf('%04d',i) '.png'])
            
        end
                    
        if SubPlot==0
           if FrameLab==0
                    %title(['Flow field motor frame, t=' num2str(time(i))])
                    if SaveSnapshot==0
                        title('Flow field motor frame')
                    end
           elseif FrameLab==1
                    %title(['Flow field lab frame, t=' num2str(time(i))])
                    %title('Flow field lab frame')
           elseif FrameLab==2
                    title(['Flow field droplet frame, t=' num2str(time(i))])
           end
        else
          title(['t=' num2str(time(i))])
        end
        
        
        
        if SaveVelField==1

            
            fig = gcf;
            fig.Position = [400 200 width height];
            name = 'BubbleFrameVelFieldMask';
            print('-dpng','-loose','-r100',[dest name sprintf('%04d',i) '.png'])
            
            
        end

        if SubPlot==1
            figure(1)
            subplot((numel(time)+2)/3,3,1)
            hold on
            plot(ttt(end),vnow(end),'.','MarkerSize',40)
        else
            figure(2)
            hold on
            plot(ttt(end),vnow(end),'.','MarkerSize',40)
        end

    end

    
        
    if numel(time)==4
            legend('v',['t=' num2str(time(1))],['t=' num2str(time(2))],['t=' num2str(time(3))],['t=' num2str(time(4))],'Location','Best')
    elseif numel(time)==3
            legend('v',['t=' num2str(time(1))],['t=' num2str(time(2))],['t=' num2str(time(3))],'Location','Best')
    end
        
    title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta) ' x_0=' num2str(PARAM.start)])

end

%plot stresses acting on the wall
if PlotStresses==1
    
    figure
    subplot((numel(time)+1)/3,3,1)
    plot(tttHere,vDrop,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('v_{drop}')
    
    for i = 1:numel(time)
        
        display([num2str(i) ' of ' num2str(numel(time))])
        
        %which iteration?
        if i==1
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1);
        else
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1)-1;
        end
        
        %elements on wall
        m = find(risa(2:end,ite)==risa(1,ite));
        indNode = find(risb(:,ite)==0,2,'first');
        q = indNode(2)-m-3;
        
        %velocity data
        vnow = vDrop(1:ite);
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(i);
        
        %wall coordinates
        aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
        
        %elements lower part of wall
        xEnd = aMotor(1)-PARAM.L*cos(PARAM.theta);
        [~,m] = min(abs(aMotor(2:end)-xEnd));
        aMotor = aMotor(1:m+1);     bMotor = bMotor(1:m+1);
        
        %stresses on wall
        FX = risy(1:2:2*m-1,ite);
        FY = risy(2:2:2*m,ite);
        
        subplot((numel(time)+1)/3,3,i+1)
        plot((aMotor(1:m)+aMotor(2:m+1))/2,FX(1:end))
        %plot((aMotor(2:m-1)+aMotor(3:m))/2,FX(2:end-1))
        grid on
        xlabel('x')
        ylabel('f_x')
        title(['t=' num2str(time(i))])
        
        %plot dot for velocity
        subplot((numel(time)+1)/3,3,1)
        hold on
        plot(ttt(end),vnow(end),'.','MarkerSize',40)
        
    end
    
end

if PlotDropVel==1
    figure
    subplot(3,1,1)
    plot(tttHere,v,'LineWidth',2)
    grid on
    xlabel('time')
    ylabel('v_{motor}')
    %hold on
    subplot(3,1,2)
    plot(tttHere,vDrop,'LineWidth',2)
    grid on
    xlabel('time')
    ylabel('v_{drop}')
    subplot(3,1,3)
    dAdt = diff(AreaDrop)./diff(tttHere);
    %plot(tttHere,AreaDrop,'LineWidth',2)
    %hold on
    plot(tttHere(1:end-1),dAdt,'LineWidth',2)
    grid on
    %title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    grid on
    xlabel('time')
    ylabel('dA/dt')
    %legend('motor velocity','bubble velocity','location','Best')
end

%zoom velocity field
if Zoom==1
    
    %which iteration?
    ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
    
    %elements on wall
    m = find(risa(2:end,ite)==risa(1,ite));
    indNode = find(risb(:,ite)==0,2,'first');
    q = indNode(2)-m-3;
    
    %wall and drop coordinates
    aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
    aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);
    
    %field points
    xSingMotor = (aMotor(1:end-1)+aMotor(2:end))/2;
    ySingMotor = (bMotor(1:end-1)+bMotor(2:end))/2;
    
    %center of mass
    xcm = center_mass(aDrop,bDrop);
    
    %select solution
    numSing = numel(aDrop)+numel(aMotor)-1;
    solution = risy(1:2*numSing+1,ite);
    
    %build grid for velocity field
    xxx = linspace(xLim(1),xLim(2),MeshZ*(xLim(2)-xLim(1)));
    yyy = linspace(yLim(1),yLim(2),MeshZ*(yLim(2)-yLim(1)));
    [X,Y] = meshgrid(xxx,yyy);
    
    %compute velocity field
    [U,V,U0] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',solution,PARAM.visc,X,Y,PARAM,FrameLab,MASKzoom,decompose);
    U = U-U0;
    
    %elements lower part of wall
    mOld = m;
    xEnd = aMotor(1)-PARAM.L*cos(PARAM.theta);
    [~,m] = min(abs(aMotor(2:end)-xEnd));
    
    %elements upper part of wall
    yFirst = bMotor(m+1)+2*PARAM.thickness*cos(PARAM.theta);
    [~,mFirst] = min(abs(bMotor-yFirst));
    mFirst = mFirst-1;
    yLast = bMotor(mFirst)+PARAM.L*sin(PARAM.theta);
    [~,mLast] = min(abs(bMotor(mFirst:end)-yLast));
    mLast = mLast-1+mFirst;
    
    %elements a portion of lower part of wall
    [~,mFirstPortion] = min(abs(aMotor(1:m+1)-xLim(2)));
    mFirstPortion = mFirstPortion-1;
    [~,mLastPortion] = min(abs(aMotor(1:m+1)-xLim(1)));
    mLastPortion = mLastPortion-1;
    
    %stresses on inner wall
    FX = risy(1:2:2*m-1,ite);
    FY = risy(2:2:2*m,ite);
    
    %stress on a portion of inner wall
    FXportion = risy(2*mFirstPortion-1:2:2*mLastPortion-1,ite);
    FYportion = risy(2*mFirstPortion:2:2*mLastPortion,ite);
    
    %compute integral of inner forces
    dl = sqrt(diff(aMotor).^2+diff(bMotor).^2);
    xInner = xSingMotor(1:m);   yInner = ySingMotor(1:m);
    dlInner = dl(1:m);
    StressInnerInt = 2*pi*sum(yInner.*dlInner.*FX)
    
    %stresses on outer wall
    FXouter = risy(2*mFirst-1:2:2*mLast-1,ite);
    FYouter = risy(2*mFirst:2:2*mLast,ite);

    %compute integral of portion forces
    xPortion = xSingMotor(mFirstPortion:mLastPortion);  yPortion = ySingMotor(mFirstPortion:mLastPortion);
    dlPortion = dl(mFirstPortion:mLastPortion);
    StressPortionInt = 2*pi*sum(yPortion.*dlPortion.*FXportion)
    StressPortionIntPercent = StressPortionInt/StressInnerInt
    
    %compute integral of portion of inner forces
    xOuter = xSingMotor(mFirst:mLast);  yOuter = ySingMotor(mFirst:mLast);
    dlOuter = dl(mFirst:mLast);
    StressOuterInt = 2*pi*sum(yOuter.*dlOuter.*FXouter)
    
    %check stress free
    StressFull = 2*pi*sum(ySingMotor.*dl.*risy(1:2:2*mOld-1,ite))
    
    figure
    %plot motor
    subplot(2,1,1)
    plot(aMotor,bMotor,'k')
    hold on
    plot(aMotor,-bMotor,'k')
    plot([aDrop; flip(aDrop)],[bDrop; -flip(bDrop)],'r')
    plot([xLim(1) xLim(2)],[yLim(1) yLim(1)],'--k')
    plot([xLim(1) xLim(1)],[yLim(1) yLim(2)],'--k')
    plot([xLim(2) xLim(2)],[yLim(1) yLim(2)],'--k')
    plot([xLim(1) xLim(2)],[yLim(2) yLim(2)],'--k')
    plot(xcm,0,'k.','MarkerSize',20)
    %quiver([Xbig Xbig],[Ybig -Ybig],[Ubig Ubig],[Vbig -Vbig])
    axis equal
    axis([-2 12 -3 3])
    xlabel('x')
    ylabel('r')
    grid on
    
    %plot velocity field
    subplot(2,1,2)
    plot(aMotor,bMotor,'k')
    hold on
    plot(aMotor,-bMotor,'k')
    plot([aDrop; flip(aDrop)],[bDrop; -flip(bDrop)],'r')
    quiver(X,Y,U,V)
    axis equal
    axis([xLim(1) xLim(2) yLim(1) yLim(2)])
    xlabel('x')
    ylabel('r')
    grid on
    
    %plot inner stresses
    figure
    subplot(2,1,1)
    plot((aMotor(1:m)+aMotor(2:m+1))/2,FX)
    grid on
    xlabel('x')
    ylabel('f_x')
    axis([xLim(1) xLim(2) min(FX)+0.5*min(FX) max(FX)+0.5*max(FX)])
    title('inner wall')
    
    %plot outer stresses
    subplot(2,1,2)
    plot(xSingMotor(mFirst:mLast),FXouter)
    grid on
    xlabel('x')
    ylabel('f_x')
    axis([xLim(1) xLim(2) min(FXouter)-0.5*min(FXouter) max(FXouter)+0.5*max(FXouter)])
    title('outer wall')
    
end

%check velocty decay
if CheckDecay==1
    
   %which iteration?
   ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
    
   %elements on wall
   m = find(risa(2:end,ite)==risa(1,ite));
   indNode = find(risb(:,ite)==0,2,'first');
   q = indNode(2)-m-3;
    
   %wall and drop coordinates
   aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
   aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);
    
   %select solution
   numSing = numel(aDrop)+numel(aMotor)-1;
   solution = risy(1:2*numSing+1,ite);
   
   %front decay
   xfront = logspace(0,dist,nDecay)+max(aMotor);
   yfront = zeros(1,numel(xfront));
   
   %compute front velocity decay
   [Ufront,Vfront,~,pFront] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',PARAM.visc,xfront,yfront,PARAM,FrameLab,MASK,decompose);
   UUUfront = sqrt(Ufront.^2+Vfront.^2);
   
   %lateral decay
   xlateral = max(aMotor)*ones(1,nDecay)-1;
   ylateral = logspace(0,dist,nDecay)+0.1;
   
   %compute lateral velocity decay
   [Ulateral,Vlateral,~,pLateral] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',PARAM.visc,xlateral,ylateral,PARAM,FrameLab,MASK,decompose);
   UUUlateral = sqrt(Ulateral.^2+Vlateral.^2);
   
   %posterior decay
   xback = -logspace(0,dist,nDecay);
   yback = zeros(1,numel(xback));
   
   %compute lateral velocity decay
   [Uback,Vback,~,pBack] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',PARAM.visc,xback,yback,PARAM,FrameLab,MASKzoom,decompose);
   UUUback = sqrt(Uback.^2+Vback.^2);
   
   %plot velocity decay
   figure
   loglog(xfront,UUUfront,'o-')
   xlabel('\rho')
   ylabel('|U|')
   grid on
   hold on
   loglog(ylateral,UUUlateral,'*-')
   loglog(abs(xback),UUUback,'x-')
   legend('front','lateral','back')
   title('Velocity decay')
   
%    figure
%    loglog(xfront,abs(pFront),'o-')
%    xlabel('\rho')
%    ylabel('|p|')
%    grid on
%    hold on
%    loglog(ylateral,abs(pLateral),'*-')
%    loglog(abs(xback),abs(pBack),'x-')
%    legend('front','lateral','back')
%    title('Pressure decay')
   
end

if ModifyStresses==1
    
    %which iteration?
    ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
    
    %elements on wall
    m = find(risa(2:end,ite)==risa(1,ite));
    indNode = find(risb(:,ite)==0,2,'first');
    q = indNode(2)-m-3;
    
    %wall and drop coordinates
    aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
    aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);
    
    %field points
    xSingMotor = (aMotor(1:end-1)+aMotor(2:end))/2;
    ySingMotor = (bMotor(1:end-1)+bMotor(2:end))/2;
    
    %stresses on wall
    FX = risy(1:2:2*m-1,ite);
    FY = risy(2:2:2*m,ite);
    
    %move everythinh with this pressure
    pResize = signPressure*min(FX);
    
    %compute normal vector
    no = sqrt(diff(aMotor).^2+diff(bMotor).^2);
    nx = -diff(bMotor)./no;  ny = diff(aMotor)./no;
    
    %resize stress by contsant pressure
    FXresize = (FX./nx+pResize).*nx;
    
%     figure
%     plot(xSingMotor,nx,xSingMotor,ny)
%     grid on
%     legend('n_x','n_y')
    
    figure
    plot(xSingMotor,FX)
    hold on
    plot(xSingMotor,FXresize,'-')
    grid on
    xlabel('x')
    ylabel('F_x')
    legend('Original','Resized','Location','Best')
    
    %COMPUTE STRESSES IN DIFFERENT LOCATION OF THE WALL
    %lower wall
    xEnd = aMotor(1)-PARAM.L*cos(PARAM.theta);
    [~,mLower] = min(abs(aMotor(2:end)-xEnd));
    
    %stresses on inner wall
    FXinner = FXresize(1:mLower);
    %FYinner = risy(2:2:2*mLower,ite);
    
    %compute integral of inner forces
    dl = sqrt(diff(aMotor).^2+diff(bMotor).^2);
    xInner = xSingMotor(1:mLower);   yInner = ySingMotor(1:mLower);
    dlInner = dl(1:mLower);
    StressInnerInt = 2*pi*sum(yInner.*dlInner.*FXinner)
    
    %elements upper part of wall
    yFirst = bMotor(mLower+1)+2*PARAM.thickness*cos(PARAM.theta);
    [~,mFirst] = min(abs(bMotor-yFirst));
    yLast = bMotor(mFirst)+PARAM.L*sin(PARAM.theta);
    [~,mLast] = min(abs(bMotor(mFirst:end)-yLast));
    mLast = mLast-2+mFirst;
    
    %compute integral of upper forces
    FXouter = FXresize(mFirst:mLast);
    xOuter = xSingMotor(mFirst:mLast);  yOuter = ySingMotor(mFirst:mLast);
    dlOuter = dl(mFirst:mLast);
    StressOuterInt = 2*pi*sum(yOuter.*dlOuter.*FXouter)
    
    %elements of left cap
    xLeft = xSingMotor(mLower+1:mFirst-1);
    yLeft = ySingMotor(mLower+1:mFirst-1);
    
    %compute integral of left forces
    FXleft = FXresize(mLower+1:mFirst-1);
    dlLeft = dl(mLower+1:mFirst-1);
    StressLeftInt = 2*pi*sum(yLeft.*dlLeft.*FXleft)
    
    %elements of right cap
    xRight = xSingMotor(mLast+1:end);
    yRight = ySingMotor(mLast+1:end);
    
    %compute integral of right forces
    FXright = FXresize(mLast+1:end);
    dlRight = dl(mLast+1:end);
    StressRightInt = 2*pi*sum(yRight.*dlRight.*FXright)
    
    %check stress free
    StressFull = 2*pi*sum(ySingMotor.*dl.*FXresize)
    
    %plot inner and outer stresses
    figure
    subplot(2,1,1)
    plot(xInner,FXinner)
    grid on
    xlabel('x')
    ylabel('f_x')
    
    subplot(2,1,2)
    plot(xOuter,FXouter)
    grid on
    xlabel('x')
    ylabel('f_x')
    
%     figure
%     plot(xInner,yInner)
%     hold on
%     plot(xOuter,yOuter)
%     plot(xLeft,yLeft)
%     plot(xRight,yRight)
%     plot(xSingMotor,ySingMotor,'o')
%     grid on
%     axis equal
%     
%     figure
%     plot(aMotor(1:mLower+1),bMotor(1:mLower+1))
%     hold on
%     plot(aMotor(mLower+1:mFirst+1),bMotor(mLower+1:mFirst+1))
%     plot(aMotor(mFirst+1:mLast+1),bMotor(mFirst+1:mLast+1))
%     plot(aMotor(mLast+1:end),bMotor(mLast+1:end))
%     grid on
%     axis equal
    
end

if ComputeDropPos==1
    
    figure
    %subplot(2,1,1)
    plot(tttHere,dropPos,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('pos_{bubble}')
    
    %subplot(2,1,2)
    figure
    plot(dropPos(1:end-1),diff(AreaDrop)./diff(dropPos),'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('space')
    ylabel('dA_{bubble}/dx')
    
end

if ComputeDropVel==1
    
    figure
    plot(tttHere,vDrop,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('V_{bubble}')
end

if ComputeDropVolume==1
    
    figure
    errV = (volumeDrop-volumeDrop(1))/volumeDrop(1);
    plot(tttHere,errV,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('errV')
end

if ComputeDropSurface==1
    
    %subplot(2,1,1)
    figure
    plot(tttHere,AreaDrop,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('A_{bubble}')
    
    if SubPlot==1
        figure(1)
        subplot(2,3,2)
    else
        figure
    end
    plot(tttHere(1:end-1),diff(AreaDrop)./diff(tttHere),'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
    hold on
    grid on
    xlabel('time')
    ylabel('dA_{bubble}/dt')
    
end

%plot pressure
if computePressureField==1 && PlotVelField
    
   figure
%    contour(X,Y,p,linspace(-0.5,0.5,1000))
%    hold on
%    contour(X,-Y,p,linspace(-0.5,0.5,1000))
   contourf(X,Y,p,[-0.5 0.5])
   hold on
   contourf(X,-Y,p)
   axis equal
   axis([xRange1 xRange2 -yRange2 yRange2])
   plot([aDrop flip(aDrop)],[bDrop -flip(bDrop)],'r-',aMotor,bMotor,'k',aMotor,-bMotor,'k','MarkerSize',2,'LineWidth',2)
   hold off
   xlabel('x')
   ylabel('r')
   title('Pressure field')
   colorbar
    
end

%plot pressure
if PlotPressureSlice==1
    
   X = linspace(xLimit(1),xLimit(2),(xLimit(2)-xLimit(1))*MeshFineSlice);
   Y = ySlice*ones(1,numel(X));
   mHalf = round(numel(X)/2);
   X0 = X(mHalf);   Y0 = Y(mHalf);
   
   %rotation
   X = (X-X0)*cos(thetaSlice) - (Y-Y0)*sin(thetaSlice) + X0;
   Y = (X-X0)*sin(thetaSlice) + (Y-Y0)*cos(thetaSlice) + Y0;
    
   %compute pressure along a line
   [U,V,~,p] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',visc,X,Y,PARAM,FrameLab,MASK,decompose);
    
   figure
   plot([aDrop flip(aDrop)],[bDrop -flip(bDrop)],'r-',aMotor,bMotor,'k',aMotor,-bMotor,'k','MarkerSize',2,'LineWidth',2)
   axis equal
   %axis([xRange1 xRange2 -yRange2 yRange2])
   hold on
   plot(X,Y)
   xlabel('x')
   ylabel('r')
   grid on
   title('Pressure along a line')
   %title('S')
   %colorbar
   
   figure
   semilogy(X,abs(p))
   xlabel('x')
   ylabel('p')
   grid on
   title('Pressure along a line')
   
   figure
   plot(X,U)
   xlabel('x')
   ylabel('u')
   grid on
   title('Axial velocity along a line')
   
   figure
   plot(X,V)
   xlabel('x')
   ylabel('v')
   grid on
   title('Radial velocity along a line')
    
end

if plotFinalCurvature==1
    
   %compute curvature
   [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(aDrop',bDrop');
   K1 = curv_spline2(bx,by,cx,cy,dx,dy);
   
   %compute normal for nodes
   N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
        -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];
    
   %compute curveture
   K2 = N(2,:)./bDrop';
   K2([1 end]) = K1([1 end]);
   K = K1+K2;
    
   figure 
   plot(aDrop,K1)
   hold on
   plot(aDrop,K2)
   plot(aDrop,K,'k-')
   grid on
   xlabel('x')
   ylabel('k')
   legend('k_1','k_2','k','Location','Best')
   
end

if checkCapillaryTwoThird==1
            
    figure
    loglog(-vDrop,hBre,'-x')
    hold on
    loglog(-vDrop,abs(vDrop).^(scaling)*mult/(abs(vDrop(1)).^(scaling)),'-k')
    %loglog(-vDrop,abs(vDrop).^(scaling),'-k')
    xlabel('Ca')
    ylabel('h')
    title('Scaling')
    grid on
    legend('Data from simuation',[num2str(scaling) ' scaling'],'Location','Best')
    
end

if plotStresses==1
    
    [~,ind] = min(abs(time(end)-(0:dt*PARAM.checkpoint:dt*PARAM.checkpoint*PARAM.loop)));
    
    %compute stresses
    y = risy(:,ind);
    fx = y(1:2:2*m-1);
    fy = y(2:2:2*m);
    
    figure
    plot((aWall(1:m)+aWall(2:m+1))/2,fx)
    hold on
    plot((aWall(1:m)+aWall(2:m+1))/2,fy)
    title('Stresses on the wall')
    grid on
    xlabel('x')
    ylabel('f_x,f_y')
    
end














