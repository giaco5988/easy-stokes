%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
PARAM.flux = 1;             % BC for flux

%shape, number of modesc
nMode = 10;
PARAM.maxCurv = 1e2;

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
plotField = 0;

%geometry parameters for Laplace
PARAM.panels = 1;                   % number of panels per block
PARAM.n = 200;                      % number of element per panel
%PARAM.typePanel = 1;               % 0 is a straight line, 1 ia an arc
PARAM.typeBClaplace = 2;            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = 0;     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = 0;     % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = 3;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 0;    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 0;    % 0 is straight, 1 is curved (spline)
PARAM.panelType = 1;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)

%function profile for BCs
PARAM.fluxBC{1} = 1;

%volume
PARAM.V0 = 4/3*r^3*pi;

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)

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

%filename
PARAM.filename = ['spherePhoretic_n=' num2str(sum(PARAM.n)) '_nMode=' num2str(nMode) '.mat'];

%finite differences
PARAM.D1{1} = finiteDifference1D(PARAM.n,[2 0],1);

%options
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolCon',1e-06,'TolX',1e-06,'TolFun',1e-06,'MaxFunEvals',1e3);

%funcion to maximize (particle speed)
fToMinimize = @(unk) computeParticleSpeed(unk,PARAM);

%funcion constraint
fConstraint = @(unk) curvatureConstraint(unk,PARAM);

%initial condition is a sphere
fModeInitial = zeros(nMode,1);

%minimize
%fMode = fminsearch(U,fModeInitial,options);

A = []; b = []; Aeq = []; beq = []; lb = -10*ones(nMode,1); ub = 10*ones(nMode,1);
fMode = fmincon(fToMinimize,fModeInitial,A,b,Aeq,beq,lb,ub,fConstraint,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot final shape
%physical grid
theta = linspace(0,pi,PARAM.n+1)';

%detrmine first mode such that the volume is conserved
fVolume = @(unk) ModifyVolumeModes2(theta,fMode,unk,PARAM.V0);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
f0 = fsolve(fVolume,1,options);

%compute shape
fMode = [f0; fMode];

figure
loglog(abs(fMode))
xlabel('n')
ylabel('f_n')
grid on

r = LegendreBuildShape(theta,fMode,0);
figure
xHere = r.*cos(theta);
yHere = r.*sin(theta);
plot(xHere,yHere,'k')
hold on
plot(xHere,-yHere,'k')
grid on
axis equal
axis([-3 3 -3 3])
xlabel('x')
ylabel('r')
title(['Shape optimization, ' num2str(numel(fMode)-1) ' modes'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solve laplace equation
x{1} = xHere';   y{1} = yHere';
conc = BEM_Laplace(x,y,PARAM);

%compute slip
vSlip = computeSlipVelPhoretic(x{1},y{1},conc(1:PARAM.n),PARAM,1,PARAM.D1{1});

%function profile for BCs
PARAM.velBC{1} = vSlip;

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);
U = yStokes(end);

display(['U0=' num2str(U)])

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









