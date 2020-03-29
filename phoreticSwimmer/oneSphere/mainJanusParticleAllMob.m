%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                          % average cone radius
PARAM.flux = [1 0];            % BC for flux
muC = 0;

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;

%geometry parameters for Laplace
nPerLenght = 50;
PARAM.panels = 2;
PARAM.n = round([acos(muC)*nPerLenght (pi-acos(muC))*nPerLenght]);                              % number of element per panel
%PARAM.typePanel = 1;                       % 0 is a straight line, 1 ia an arc
PARAM.typeBClaplace = [2 2];                % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0];         % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0];         % 0 is straight, 1 is curved (spline)

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                  % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;
PARAM.deflationBlock = 1;

%function profile for BCs
PARAM.fluxBC{1} = PARAM.flux(1);
PARAM.fluxBC{2} = PARAM.flux(2);

%build geometry
[x{1},y{1}] = buildArc(r,0,0,[0 acos(muC)],PARAM,1);
[x{2},y{2}] = buildArc(r,0,0,[acos(muC) pi],PARAM,2);

%plot conical motor shape
% figure
% plot(x{1},y{1},'r')
% hold on
% plot(x{2},y{2},'k')
% plot(x{2},-y{2},'k')
% plot(x{1},-y{1},'r')
% grid on
% xlabel('x')
% ylabel('r')
% axis equal
% title('Geometry')

%print to screen
printToScreenLaplace(PARAM)
printToScreenStokes(PARAM)

%filename
PARAM.filename = ['janusSphere_n=' num2str(sum(PARAM.n)) '.mat'];

%finite differences
D1{1} = finiteDifference1D(PARAM.n(1),[2 0],1);
D1{2} = finiteDifference1D(PARAM.n(2),[2 0],1);

%solve laplace equation
conc = BEM_Laplace(x,y,PARAM);

%compute slip
[vSlip1,l1] = computeSlipVelPhoretic(x{1},y{1},conc(1:PARAM.n(1)),PARAM,1,D1{1});
[vSlip2,l2] = computeSlipVelPhoretic(x{2},y{2},conc(PARAM.n(1)+1:end),PARAM,2,D1{2});

%function profile for BCs
PARAM.velBC{1} = vSlip1;
PARAM.velBC{2} = vSlip2;

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%plot stresses
figure
fx = yStokes(1:2:end-2);
fy = yStokes(2:2:end-1);
plot([l1; l1(end)+l2],fx)
hold on
plot([l1; l1(end)+l2],fy)
grid on
xlabel('x')
ylabel('f_x,f_y')
title('Stresses on wall')

%velocity
U0 = yStokes(end);
Uan = 1/4*(1-muC^2)*PARAM.flux(1);
err = abs(U0-Uan)./abs(Uan);

display(['U=' num2str(U0)])
display(['errU=' num2str(err)])

%stresslet
%S = computeStressLetAxis(x,y,[vSlip1; vSlip2],1:2,PARAM);
S = 5*pi*muC*(1-muC^2);
Snum = computeStressLetAxisNumerical(x,y,yStokes(1:end-1),PARAM);

display(['S=' num2str(S)])
display(['Snum=' num2str(Snum)])

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









