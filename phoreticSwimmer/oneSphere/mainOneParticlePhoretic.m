%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
PARAM.flux = -1;             % BC for flux

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
PARAM.deflation = 0;

%geometry parameters for Laplace
PARAM.panels = 1;
PARAM.n = 100;                              % number of element per panel
%PARAM.typePanel = 1;                        % 0 is a straight line, 1 ia an arc
PARAM.typeBClaplace = 2;            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = 0;     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = 0;     % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = 3;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 0;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 0;    % 0 is straight, 1 is curved (spline)
PARAM.panelType = 1;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;
PARAM.deflationBlock = 1;

%function profile for BCs
theta = linspace(0,pi,PARAM.n+1)';
theta = (theta(1:end-1)+theta(2:end))/2;
flux = (cos(theta)+1)/2;
PARAM.fluxBC{1} = flux;

%build geometry
[x{1},y{1}] = buildArc(r,0,0,[0 pi],PARAM,1);

%plot conical motor shape
figure
plot(x{1},y{1},'k')
hold on
plot(x{1},-y{1},'k')
grid on
xlabel('x')
ylabel('r')
axis equal
title('Geometry')

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)

%filename
PARAM.filename = ['spherePhoretic_n=' num2str(sum(PARAM.n)) '.mat'];

%finite differences
D1{1} = finiteDifference1D(PARAM.n,[2 0],1);

%solve laplace equation
conc = BEM_Laplace(x,y,PARAM);

%compute slip
[vSlip,l] = computeSlipVelPhoretic(x{1},y{1},conc(1:PARAM.n),PARAM,1,D1{1});

%function profile for BCs
PARAM.velBC{1} = vSlip;

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%plot stresses
figure
fx = yStokes(1:2:end-2);
fy = yStokes(2:2:end-1);
plot(l,fx)
hold on
plot(l,fy)
grid on
xlabel('x')
ylabel('f_x,f_y')
title('Stresses on wall')

display(['U0=' num2str(yStokes(end))])

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









