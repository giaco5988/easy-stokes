%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
%close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                          % average cone radius
rBubble = 0.022;                  % radius of the bubble
L = 10;                         % motor lenght
h = 0.2;                        % wall thickness
thetaCone = 0;                      % inclinitation of the cone
InPos = 0;                   % initial bubble position
K = 0.1;
beta = 1;
%Ad3 = 1;

%time marching
loop = 1000;
dt = 0.001;

rBubble0 = rBubble;

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
plotField = 0;
threeD = 1; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.5];

%geometry parameters
nPerLenght = 10;
PARAM.panels = [4 1];                              % panels per block
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) round(10*rBubble*pi/2*nPerLenght)];                  % number of element per panel
PARAM.geometryPanel = [0 1 0 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2 rBubble];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L InPos];
PARAM.y0_Circle = [nan r nan r 0];
PARAM.SPlinesType = [nan nan nan nan nan];

%function profile for BCs
PARAM.fluxBC{1} = 0;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = -1;
PARAM.fluxBC{4} = 0;
PARAM.fluxBC{5} = 0;
PARAM.concBC{1} = 0;
PARAM.concBC{2} = 0;
PARAM.concBC{3} = 0;
PARAM.concBC{4} = 0;
PARAM.concBC{5} = (2/rBubble + beta)*K;

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [0 0 0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)

%Laplace
PARAM.typeBClaplace = [2 2 2 2 1];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0];     % 0 is straight, 1 is curved (spline)

%build geometry
xC = 0;   yC = r;
tParametric = parametricGrid(PARAM);

%build physical shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

%rigid body rotation
[x,y,tParametric,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd] = rigidBodyRotation(x,y,tParametric,xC,yC,thetaCone,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd,PARAM.geometryPanel,1:4);

%plot conical motor shape
figure
plotGeometry(x,y,threeD,theta3D,color3D,transp,PARAM)
hold off
grid on
xlabel('x')
ylabel('r')
axis equal
title('Motor geometry')

%print to screen
printToScreenLaplace(PARAM)
%printToScreenStokes(PARAM)

%filename
PARAM.filename = ['coneExpandBubble_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(thetaCone) '_h=' num2str(h) '_rBubble=' num2str(rBubble) '_InPos=' num2str(L/2) '.mat'];

%finite differences
D1{1} = finiteDifference1D(PARAM.n(3),[2 0],1);

for i = 1:loop
    
    display(['T=' num2str(dt*(i-1))])

    %solve laplace equation
    [yLaplace,Xsing,Ysing] = BEM_Laplace(x,y,PARAM);

    %compute slip
    [vSlip,l] = computeSlipVelPhoretic(x{3},y{3},yLaplace(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{1});

    %compute flow rate
    Q = computeFlowRate(x{5},y{5},yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n)),5,PARAM);
    manyQ(i) = Q;
    manyV(i) =  axis_int_gauss_vect(x{5},y{5});
    
    %inflate bubble
    Un = Q/(3*beta^rBubble^2+4*rBubble);
    rBubble = rBubble + dt*Un;
    manyR(i) = rBubble;
    
    %new concentration on the bubble
    PARAM.concBC{5} = K*(beta + 2/rBubble);

    %function profile for BCs
    PARAM.velBC{1} = 0;
    PARAM.velBC{2} = 0;
    PARAM.velBC{3} = vSlip;
    PARAM.velBC{4} = 0;

    %solve Stokes equation
    %[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);
    
    %new bubble geometry
    %PARAM.n(5) = round(10*rBubble*pi/2*nPerLenght);
    [x{5},y{5}] = buildArc(rBubble,InPos,0,[0 pi],PARAM,5);
    
    %bubble disapperas
    if min(y{5}<0.009)
        display('Bubble has vanished')
        break;
    end
    
    %bubble touches than the tube
    phi = FigInOut(x{5},y{5},[x{1} x{2} x{3} x{4}]',[y{1} y{2} y{3} y{4}]');
    comp = max(phi>pi);
    if comp
        display('Bubble touches the tube')
        break;
    end
    
    %plot new geometry
    figure(1)
    plotGeometry(x,y,threeD,theta3D,color3D,transp,PARAM)
    hold off
    grid on
    xlabel('x')
    ylabel('r')
    axis equal
    %axis(10*[-0.5 1.5 -0.5 0.5])
    title('Motor geometry')
    axis off
    drawnow

end

%display(['U0=' num2str(yStokes(end))])

display('The end')

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

figure
plot(0:dt:(i-2)*dt,manyQ(1:end-1))
%semilogy(0:dt:(i-1)*dt,manyQ-manyQ(1))
xlabel('t')
ylabel('Q')
grid on
title(['Oxygen flux z_0=' num2str(InPos) ' r_0=' num2str(rBubble0)])

figure(3)
hold on
%plot(0:dt:(i-2)*dt,manyR(1:end-1),'-')
semilogy(0:dt:(i-2)*dt,manyR(1:end-1),'-')
%semilogy(0:dt:(i-1)*dt,manyV-manyV(end),'-')
xlabel('t')
ylabel('r')
grid on
title(['Bubble radius z_0=' num2str(InPos) ' r_0=' num2str(rBubble0)])








