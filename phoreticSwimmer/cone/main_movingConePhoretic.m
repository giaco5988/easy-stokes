%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 0.1;                      % average cone radius
L = 1;                   % motor lenght
h = 0.05;                   % wall thickness
theta = pi/16;               % inclinitation of the cone

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
PARAM.deflationBlock = 1;

%geometry parameters
nPerLenght = 500;
PARAM.panels = 4;
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel
PARAM.typePanel = [0 1 0 1];                % 0 is a straight line, 1 ia an arc
PARAM.typeBClaplace = [2 2 2 2];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0];     % 0 is straight, 1 is curved (spline)

%function profile for BCs
PARAM.fluxBC{1} = 0;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = -1;
PARAM.fluxBC{4} = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [1 1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;

%build geometry
xC = 0;   yC = r;
[x1,y1] = buildStraightLine(L,r+h/2,0,r+h/2,PARAM,1);
[x2,y2] = buildArc(h/2,0,r,[pi/2 3*pi/2],PARAM,2);
[x3,y3] = buildStraightLine(0,r-h/2,L,r-h/2,PARAM,3);
[x4,y4] = buildArc(h/2,L,r,[-pi/2 pi/2],PARAM,4);
x{1} = x1;  x{2} = x2;  x{3} = x3;  x{4} = x4;
y{1} = y1;  y{2} = y2;  y{3} = y3;  y{4} = y4;

%solid body rotation
[x,y] = rigidBodyRotation(x,y,xC,yC,theta,1:4);

%plot conical motor shape
plotGeometry(x,y,PARAM)
grid on
xlabel('x')
ylabel('r')
axis equal
title('Motor geometry')

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)

%filename
PARAM.filename = ['conePhoretic_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

%finite differences
D1{1} = finiteDifference1D(PARAM.n(3),[2 0],1);

%solve laplace equation
yLaplace = BEM_Laplace(x,y,PARAM);

%compute slip
[vSlip,l] = computeSlipVelPhoretic(x{3},y{3},yLaplace(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{1});

%plot concentration
figure
plot(l,yLaplace(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),'x-')
hold on
plot(l,vSlip,'x-')
xlabel('l')
ylabel('c,\nabla c')
grid on
title('Surface concentration and derivative')

%function profile for BCs
PARAM.velBC{1} = 0;
PARAM.velBC{2} = 0;
PARAM.velBC{3} = vSlip;
PARAM.velBC{4} = 0;

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%compute stresslet numerically
Snum = computeStressLetAxisNumerical(x,y,yStokes(1:end-1),PARAM);

display(['U0=' num2str(yStokes(end))])
display(['S=' num2str(Snum)])

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









