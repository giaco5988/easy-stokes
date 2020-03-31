%plot velocity field

clear variables
close all

PARAM.res = '~/Documents/MATLAB/droplet_simulations/results';
%PARAM.res = '~/Documents/MATLAB/droplet_simulations/results/edgeTracking/forPaper/data';

%options
plotShape = 1;
plotVelField = 1;   streamSLICE = 1;
plotPressureField = 0;
plotCurvature = 0;
plotPressureAxis = 1;

%parametes
PARAM.Ca = -0.2;
PARAM.visc = 1;
PARAM.algorithm = 1;
V0 = 4/3*pi;
PARAM.placeShapeXCM = 1;
PARAM.remeshStart = 0;
PARAM.BC = 1;
PARAM.dropFrame = 0;
PARAM.dh = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Grid = 200;
X0 = linspace(-2,2,Grid);
Y0 = linspace(0,2,Grid);
%y0 = 2;
[X0,Y0] = meshgrid(X0,Y0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% SPACE DISCRETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PARAM.legendre = 1;                                         % choose 0 Chebyshev, 1 Legendre, 2 is Legendre-Lobatto
PARAM.dealiasing = 100;  ratioGrid = 2;                     % dealiasing
PARAM.tresholdCoeffs = 2e-3;     PARAM.Rescale = 1;         % option for accuracy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PARAM.dealiasing==98
    UnstOrSt = 'Stable';
elseif PARAM.dealiasing==100
    UnstOrSt = 'Unstable';
else
        error('No file')
end
%upload = load(['/Users/Giacomo/Documents/MATLAB/Postprocessing/oneDroplet/curvBaseState/curv' UnstOrSt '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '.mat']);
%Kbase = upload.K;
%upload = load(['/Users/Giacomo/Documents/MATLAB/Postprocessing/oneDroplet/curvBaseState/pressure' UnstOrSt '_Ca=' num2str(PARAM.Ca) '_visc=' num2str(PARAM.visc) '.mat']);
%pressureAxisBase = upload.pressureAxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose where to get the shape
PARAM.uploadShape = 2;                                       % 0 the initial shape is analytical, 1 is from continuation, 2 the initial shape is from newton Method,  3 is from edge tracking, 4 is from DNS, 5 is from Minimal Seed
PARAM.D = 0;                                              % deformation parameter of the initial condition

% in case the initial shape is analytical, choose how to draw it
PARAM.shapeEllipse = 1;                                      % 0 initial shape is sphere plus P2 Legendre, 1 initial shape is an ellipse, 2 is pertubed with 2 legendre polynomials P2 P4, 3 is asymmetric shape by P3 Legendre, 4 is P2, P3
PARAM.f2 = 1.5;                                              % if I use first two leegndre polynomials

% if the shape is uploaded from previous simulation
PARAM.elemUpload = PARAM.dealiasing;                                       % number of elemnt of uploaded data
PARAM.BCupload = 1;                                          % BC of the uploaded shape
PARAM.legendreUpload = 1;                                    % spectral basis of the previous simulation
PARAM.overlapMode = 2;  PARAM.whichmode = 2;                    % 1 perturb with legendre polynomial, 2 perturb shape with eigenmode

% if it start from continuation, newton or edge tracking
PARAM.CaUpUpload = 0.13; PARAM.CaDownUpload = -1;              % for uploading the file
PARAM.Dupload = 2.2;                                           % choose defromation parameter on bifurcation diagram
PARAM.CaUpload = PARAM.Ca;                                     % choose Ca on bifurcation diagram
PARAM.viscUpload = PARAM.visc;                                  % choose viscosity ratio for upload

% if the shape is uploaded from edge tracking or DNS
PARAM.edgeLoopUpload = 200;                                   % loops of the edge tracking simulation
PARAM.deltaEdgeUpload = 1e-6;                                % delta edge of the uloaded data
PARAM.dtUpload = 2e-3;                                       % time step of the uploaded simulation
PARAM.volCorrUpload = 0;                                     % if the simulation was using volume correction
PARAM.whichLoopUpload = 20;                                  % which loop of the edge tracking has to be loaded
PARAM.Tupload = 50;                                         % choose at which time upload shape

% if shape is uploaded from DNS
PARAM.TendUpload = 70;                                         % choose at which time upload shape

% if shape is uploaded from Minimal Seed
PARAM.A0upload = 0.2;                                         % energy amplification of previous optimization
PARAM.ThorizonUpload = 10;                                   % time horizon of previous optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of grid points
PARAM.n = round(ratioGrid*PARAM.dealiasing)+1;
PARAM.V0 = V0;
%PARAM.Tend = Tend;
%if Tstart~=PARAM.Tupload && PARAM.uploadShape==4
 %   error('Upload and start time have to be the same');
%end

%compute integration weight and differentiation matrices
if PARAM.legendre==1
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreIntDiff(0,1,PARAM.n+1);
    [PARAM.tcheb,PARAM.D1cheb,PARAM.D2cheb,PARAM.WGcheb] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==0
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.TTT] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==2
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreLobattoIntDiff(0,1,PARAM.n+1);
end

%upload shape
here = pwd;
cd('/Users/Giacomo/Documents/MATLAB/droplet_simulations/dropSpectral')
[~,x,y,~,~,~,~,PARAMupload] = initialConditionDrop(PARAM);
cd(here)

%compute velocity field
if plotVelField==1 || plotPressureField==1
    [u,v,pressure,K] = computeVelPressFieldOneDropletSpectral(x,y,X0,Y0,PARAMupload);
end

%compute velocity field along the axis
if PARAM.dealiasing==100 || PARAM.dealiasing==50
    xAxis = linspace(-1.8,1.8,100);
elseif PARAM.dealiasing==98 || PARAM.dealiasing==48
    xAxis = linspace(-1,1,100);
end
yAxis = zeros(1,numel(xAxis));
[uAxis,vAxis,pressureAxis] = computeVelPressFieldOneDropletSpectral(x,y,xAxis,yAxis,PARAMupload);

%plot vel field
if plotVelField==1
    figure
    [~,h1] = contourf(X0,Y0,sqrt(u.^2+v.^2),500);
    colorbar
    hold on
    [~,h2] = contourf(X0,-Y0,sqrt(u.^2+v.^2),500);
    set(h1,'LineColor','none')
    set(h2,'LineColor','none')
    plot(x,y,'k')
    hold on
    plot(x,-y,'k')
    title('Velocity field')
    axis equal
    xlabel('z')
    ylabel('r')
    axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])
    drawnow
    
    if streamSLICE==1

    hStream1 = streamslice(X0,Y0,u,v,0.5);
    hStream2 = streamslice(X0,-Y0,u,-v,0.5);

    set(hStream1,'LineWidth',1.5,'Color','white')
    set(hStream2,'LineWidth',1.5,'Color','white')

    elseif streamSLICE==0

        if PARAM.dealiasing==98

            %Close look, unconfined bubble
            startX1 = linspace(-2.5,2.5,10);    startX1 = startX1(2:end-1);  startY1 = 1*ones(1,numel(startX1));
            startX2 = [-0.5 0.5];  startY2 = 0.5*ones(1,numel(startX2));
            startX3 = [-0.2 0.2];  startY3 = 0.5*ones(1,numel(startX3));
            %startX3 = [];   startY3 = [];

            optionStream2 = [0.1 1000];
            optionStream3 = [0.1 2300];
            hStream1 = streamline(X0,Y0,u,v,startX1,startY1);
            hStream1bis = streamline(X0,Y0,u,v,startX2,startY2,optionStream2);
            hStream1tris = streamline(X0,Y0,u,v,startX3,startY3,optionStream3);
            hStream2 = streamline(X0,-Y0,u,-v,startX1,-startY1);
            hStream2bis = streamline(X0,-Y0,u,-v,startX2,-startY2,optionStream2);
            hStream2tris = streamline(X0,-Y0,u,-v,startX3,-startY3,optionStream3);
            xyCoord = stream2(X0,Y0,u,v,[startX1 startX2 startX3],[startY1 startY2 startY3]);

            set(hStream1,'LineWidth',1.5,'Color','white')
            set(hStream2,'LineWidth',1.5,'Color','white')
            set(hStream1bis,'LineWidth',1.5,'Color','white')
            set(hStream2bis,'LineWidth',1.5,'Color','white')
            set(hStream1tris,'LineWidth',1.5,'Color','white')
            set(hStream2tris,'LineWidth',1.5,'Color','white')

        elseif PARAM.dealiasing==100

            %Close look, unconfined bubble
            startX1 = linspace(-2.5,2.5,10);    startX1 = startX1(2:end-1);  startY1 = 1*ones(1,numel(startX1));
            startX2 = [-1 1];  startY2 = 0.5*ones(1,numel(startX2));
            startX3 = [-0.5 0.5];  startY3 = 0.5*ones(1,numel(startX3));
            %startX3 = [];   startY3 = [];

            optionStream2 = [0.1 2000];
            optionStream3 = [0.1 2300];
            hStream1 = streamline(X0,Y0,u,v,startX1,startY1);
            hStream1bis = streamline(X0,Y0,u,v,startX2,startY2,optionStream2);
            hStream1tris = streamline(X0,Y0,u,v,startX3,startY3,optionStream3);
            hStream2 = streamline(X0,-Y0,u,-v,startX1,-startY1);
            hStream2bis = streamline(X0,-Y0,u,-v,startX2,-startY2,optionStream2);
            hStream2tris = streamline(X0,-Y0,u,-v,startX3,-startY3,optionStream3);
            xyCoord = stream2(X0,Y0,u,v,[startX1 startX2 startX3],[startY1 startY2 startY3]);

            set(hStream1,'LineWidth',1.5,'Color','white')
            set(hStream2,'LineWidth',1.5,'Color','white')
            set(hStream1bis,'LineWidth',1.5,'Color','white')
            set(hStream2bis,'LineWidth',1.5,'Color','white')
            set(hStream1tris,'LineWidth',1.5,'Color','white')
            set(hStream2tris,'LineWidth',1.5,'Color','white')

        end
        
        %add arrows to stremlines
            for zzz = 1:numel(xyCoord)
                
                firstArrow = 500;
                       
                      xyCoordHere = xyCoord{zzz};
                      %xCoord = xyCoordHere(50:240:end,1);
                      %yCoord = xyCoordHere(50:240:end,2);
                      sizeHere = size(xyCoordHere);
                      xCoord = xyCoordHere(ceil(2*sizeHere(1)/3),1);
                      yCoord = xyCoordHere(ceil(2*sizeHere(1)/3),2);
                      
                      %compute velocity
                      [uLine,vLine,pressureLine] = computeVelPressFieldOneDropletSpectral(x,y,xCoord,yCoord,PARAMupload);
                      
                      dx = uLine./sqrt(uLine.^2+vLine.^2)/100;
                      dy = vLine./sqrt(uLine.^2+vLine.^2)/100;
                      dlHere = sqrt(diff(xCoord).^2+diff(yCoord).^2);
                      distOrigin = sqrt(xCoord.^2+yCoord.^2);
                      for sss = 1:numel(xCoord)
                          
                          if sss==1 && distOrigin(sss)>0
                              arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','w')
                              arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','w')
                          end
                          
                          if sss>1
                              if dlHere(sss-1)>firstArrow && distOrigin(sss)>0
                                arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','r')
                                arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','r')
                              end
                          end
                      
                      end
                      
            end
    
    end
                       
end

%plot pressure field
if plotPressureField==1
    figure
    [~,h1] = contourf(X0,Y0,pressure,500);
    colorbar
    hold on
    [~,h2] = contourf(X0,-Y0,pressure,500);
    set(h1,'LineColor','none')
    set(h2,'LineColor','none')
    plot(x,y,'k')
    hold on
    plot(x,-y,'k')
    title('Pressure field')
    axis equal
    xlabel('z')
    ylabel('r')
    axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])
    if PARAM.dealiasing==100
        caxis([0 4])
        colorbar('Ticks',[0 2 4])
    elseif PARAM.dealiasing==98
        caxis([0 3])
        colorbar('Ticks',[0 1.5 3])
    end
    colormap(inferno())
    drawnow
    
    if streamSLICE==0

        if PARAM.dealiasing==98

            %Close look, unconfined bubble
            startX1 = linspace(-2.5,2.5,10);    startX1 = startX1(2:end-1);  startY1 = 1*ones(1,numel(startX1));
            startX2 = [-0.5 0.5];  startY2 = 0.5*ones(1,numel(startX2));
            startX3 = [-0.2 0.2];  startY3 = 0.5*ones(1,numel(startX3));
            %startX3 = [];   startY3 = [];

            optionStream2 = [0.1 1000];
            optionStream3 = [0.1 2300];
            hStream1 = streamline(X0,Y0,u,v,startX1,startY1);
            hStream1bis = streamline(X0,Y0,u,v,startX2,startY2,optionStream2);
            hStream1tris = streamline(X0,Y0,u,v,startX3,startY3,optionStream3);
            hStream2 = streamline(X0,-Y0,u,-v,startX1,-startY1);
            hStream2bis = streamline(X0,-Y0,u,-v,startX2,-startY2,optionStream2);
            hStream2tris = streamline(X0,-Y0,u,-v,startX3,-startY3,optionStream3);
            xyCoord = stream2(X0,Y0,u,v,[startX1 startX2 startX3],[startY1 startY2 startY3]);

            set(hStream1,'LineWidth',1.5,'Color','white')
            set(hStream2,'LineWidth',1.5,'Color','white')
            set(hStream1bis,'LineWidth',1.5,'Color','white')
            set(hStream2bis,'LineWidth',1.5,'Color','white')
            set(hStream1tris,'LineWidth',1.5,'Color','white')
            set(hStream2tris,'LineWidth',1.5,'Color','white')

        elseif PARAM.dealiasing==100

            %Close look, unconfined bubble
            startX1 = linspace(-2.5,2.5,10);    startX1 = startX1(2:end-1);  startY1 = 1*ones(1,numel(startX1));
            startX2 = [-1 1];  startY2 = 0.5*ones(1,numel(startX2));
            startX3 = [-0.5 0.5];  startY3 = 0.5*ones(1,numel(startX3));
            %startX3 = [];   startY3 = [];

            optionStream2 = [0.1 2000];
            optionStream3 = [0.1 2300];
            hStream1 = streamline(X0,Y0,u,v,startX1,startY1);
            hStream1bis = streamline(X0,Y0,u,v,startX2,startY2,optionStream2);
            hStream1tris = streamline(X0,Y0,u,v,startX3,startY3,optionStream3);
            hStream2 = streamline(X0,-Y0,u,-v,startX1,-startY1);
            hStream2bis = streamline(X0,-Y0,u,-v,startX2,-startY2,optionStream2);
            hStream2tris = streamline(X0,-Y0,u,-v,startX3,-startY3,optionStream3);
            xyCoord = stream2(X0,Y0,u,v,[startX1 startX2 startX3],[startY1 startY2 startY3]);

            set(hStream1,'LineWidth',1.5,'Color','white')
            set(hStream2,'LineWidth',1.5,'Color','white')
            set(hStream1bis,'LineWidth',1.5,'Color','white')
            set(hStream2bis,'LineWidth',1.5,'Color','white')
            set(hStream1tris,'LineWidth',1.5,'Color','white')
            set(hStream2tris,'LineWidth',1.5,'Color','white')

        end
        
        %add arrows to stremlines
            for zzz = 1:numel(xyCoord)
                
                firstArrow = 500;
                       
                      xyCoordHere = xyCoord{zzz};
                      %xCoord = xyCoordHere(50:240:end,1);
                      %yCoord = xyCoordHere(50:240:end,2);
                      sizeHere = size(xyCoordHere);
                      xCoord = xyCoordHere(ceil(2*sizeHere(1)/3),1);
                      yCoord = xyCoordHere(ceil(2*sizeHere(1)/3),2);
                      
                      %compute velocity
                      [uLine,vLine,pressureLine] = computeVelPressFieldOneDropletSpectral(x,y,xCoord,yCoord,PARAMupload);
                      
                      dx = uLine./sqrt(uLine.^2+vLine.^2)/100;
                      dy = vLine./sqrt(uLine.^2+vLine.^2)/100;
                      dlHere = sqrt(diff(xCoord).^2+diff(yCoord).^2);
                      distOrigin = sqrt(xCoord.^2+yCoord.^2);
                      for sss = 1:numel(xCoord)
                          
                          if sss==1 && distOrigin(sss)>0
                              arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','w')
                              arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','w')
                          end
                          
                          if sss>1
                              if dlHere(sss-1)>firstArrow && distOrigin(sss)>0
                                arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','r')
                                arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','r')
                              end
                          end
                      
                      end
                      
            end
    
    end
    
end

if plotShape==1
    figure(4)
    plot(x,y,'k')
    hold on
    plot(x,-y,'k')
    axis equal
    xlabel('z')
    ylabel('r')
    axis([-3 3 -2 2])
    drawnow
end

if plotCurvature==1

    %plot curvature
    figure(5)
    plot(x,K)
    xlabel('z')
    ylabel('K')
    grid on
    title('Curvature')

    %plot curvature minus base state
    figure(6)
    plot(x,K-Kbase)
    xlabel('z')
    ylabel('\Delta K')
    grid on
    title('Excess Curvature')

end

if plotPressureAxis==1

    %plot curvature
    figure(7)
    plot(xAxis,pressureAxis)
    xlabel('z')
    ylabel('p')
    grid on
    title('Pressure on the axis')

    %plot curvature minus base state
    figure(8)
    plot(xAxis,pressureAxis-pressureAxisBase)
    %plot(xAxis,(pressureAxis-pressureAxisBase)./pressureAxisBase)
    xlabel('z')
    ylabel('\Delta p')
    grid on
    title('Excess Pressure')

end




























