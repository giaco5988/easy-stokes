%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
rBubble = 1.1;              % radius of the bubble
rBubble2 = 1;              % radius of the bubble
h = 0.2;                    % wall thickness
temp2 = 2/180*pi;    % inclinitation of the cone  
temp1 = 10;  
xcmBubble = 9.9;
xcmBubble2 = 8;

for k = 1:numel(temp2)
for i = 1:numel(temp1)
    
thetaCone = temp2(k);
L = temp1(i);

%nondimensional number
Hcc = 0.1;
beta = 10;

%rBubble = 2*Hcc/(35-Hcc);

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
plotField = 0;
threeD = 1; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.3];

%geometry parameters
nPerLenght = 50;
PARAM.panels = [4 1];                              % panels per block
%PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) round(rBubble*pi/2*nPerLenght) round(rBubble2*pi/2*nPerLenght)];                  % number of element per panel
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) round(rBubble*pi/2*nPerLenght)+20];                  % number of element per panel
%PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel
PARAM.typePanel = [0 1 0 1 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.geometryPanel = [0 1 0 1 1 1];
PARAM.xStart = [L nan 0 nan nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi pi];               % theta starting point for the arc
PARAM.rArc = [nan h/2 nan h/2 rBubble rBubble2];               % theta starting point for the arc
PARAM.x0_Circle = [nan 0 nan L xcmBubble xcmBubble2];
PARAM.y0_Circle = [nan r nan r 0 0];

%function profile for BCs
PARAM.flux = [0 0 -1 0 0];   % BC for flux
PARAM.fluxBC{1} = 0;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = PARAM.flux(3);
PARAM.fluxBC{4} = 0;
PARAM.fluxBC{5} = 0;
PARAM.fluxBC{6} = 0;
PARAM.concBC{1} = 0;
PARAM.concBC{2} = 0;
PARAM.concBC{3} = 0;
PARAM.concBC{4} = 0;
PARAM.concBC{5} = Hcc*(beta+2/rBubble);
PARAM.concBC{6} = Hcc*(beta+2/rBubble2);

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [0 0 0 0 0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.SPlinesType = [0 0 0 0 0 0];            % 1 is natural splines, 2 is symmetric on the axis

%numerics Laplace
PARAM.typeBClaplace = [2 2 2 2 1 1];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0 0];     % 0 is straight, 1 is curved (spline)

%build geometry
xC = 0;   yC = r;
tParametric = parametricGrid(PARAM);

%build physical shape
[x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

%rigid body rotation
[x,y,tParametric,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd] = rigidBodyRotation(x,y,tParametric,xC,yC,thetaCone,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd,PARAM.geometryPanel,1:4);

%base shape, parametric grid
tParametricBase = tParametric;

%plot conical motor shape
if i==1 && k==1
    figure
    plotGeometry(x,y,threeD,theta3D,color3D,transp,PARAM);
    grid on
    xlabel('x')
    ylabel('r')
    axis equal
    title('Motor geometry')
end

%print to screen
printToScreenLaplace(PARAM)
%printToScreenStokes(PARAM)

%filename
PARAM.filename = ['coneWithBubble_n=' num2str(nPerLenght) '_nBubble=' num2str(numel(PARAM.n)-4) '_xi=' num2str(L) '_theta=' num2str(thetaCone) '_h=' num2str(h) '_rBubble=' num2str(rBubble) '_InPos=' num2str(xcmBubble) '.mat'];

%finite differences
D1{1} = finiteDifference1D(PARAM.n(3),[2 0],1);

%solve laplace equation
[conc,Xsing,Ysing] = BEM_Laplace(x,y,PARAM);

if numel(PARAM.n)==5
    Q = computeFlowRate(x{5},y{5},conc(sum(PARAM.n(1:4))+1:sum(PARAM.n)),5,PARAM);
    display(['Q=' num2str(Q)])
elseif numel(PARAM.n)==6
    Q1 = computeFlowRate(x{5},y{5},conc(sum(PARAM.n(1:4))+1:sum(PARAM.n(1:5))),5,PARAM);
    Q2 = computeFlowRate(x{6},y{6},conc(sum(PARAM.n(1:5))+1:sum(PARAM.n)),6,PARAM);
    display(['Q1=' num2str(Q1) ' and Q2=' num2str(Q2)])
end

%compute slip
%[vSlip,l] = computeSlipVelPhoretic(x{3},y{3},conc(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{1});

%function profile for BCs
% PARAM.velBC{1} = 0;
% PARAM.velBC{2} = 0;
% PARAM.velBC{3} = vSlip;
% PARAM.velBC{4} = 0;

%solve Stokes equation
%[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%display(['U0=' num2str(yStokes(end))])

display('The end')

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

end
end









