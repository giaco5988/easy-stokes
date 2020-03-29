%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%physical parameters
x0 = 0;                 % center of the particle
r = 1;                  % radius of the particle
PARAM.phi = [0 0];      % this is not used
PARAM.flux = 1;     % BC for flux

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;

%geometry parameters
PARAM.panels = 1;                % number of panels per block
PARAM.n = 40;                    % number of element per panel
PARAM.typePanel = 1;         % 0 is a straight line

%numerics parameters for Laplace
PARAM.typeBClaplace = [2 2];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0];     % 0 is straight, 1 is curved (spline)

%functions for BCs
theta = linspace(0,pi,PARAM.n+1)';
theta = (theta(1:end-1)+theta(2:end))/2;
flux = (cos(theta)+1)/2;
PARAM.fluxBC{1} = flux;

%finite differences
D1{1} = finiteDifference1D(PARAM.n(1),[2 0],1);

%build geometry
[x1,y1] = buildArc(r,x0,0,[0 pi],PARAM,1);
x{1} = x1;  %x{2} = x2;
y{1} = y1;  %y{2} = y2;

printToScreenLaplace(PARAM);

%solve laplace equation
[yOut,Xsing] = BEM_Laplace(x,y,PARAM);

%compute slip
[vSlip,l] = computeSlipVelPhoretic(x{1},y{1},yOut(1:PARAM.n(1)),PARAM,1,D1{1});

%check slip velocity
[cAn,vSlipAn] = plotSpeciesFieldAnalyticalFunSmooth2(1,theta);

figure
plot(Xsing,vSlipAn)
hold on
grid on
plot(Xsing,vSlip,'--')
xlabel('\mu')
legend('Analytical','Numerical','Location','Best')
ylabel('\nabla c')

figure
err = abs(vSlip-vSlipAn)./abs(vSlipAn);
semilogy(Xsing,err)
xlabel('\mu')
ylabel('err')
grid on
title(['error n=' num2str(PARAM.n)])

%numerics parameters for Stokes
PARAM.typeBCstokes = 3;                                   % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = 1;                            % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = 0;                            % 0 is straight, 1 is curved (spline)
PARAM.velBC{1} = [vSlip; vSlip(end)];
%PARAM.velBC{1} = vSlip;
PARAM.stressBC{1} = 0;%      PARAM.stressBC{2} = 0;
PARAM.panelType = 1;                                % 0 is fixed wall, 1 is moving wall, 2 is droplet (this decide the coefficient of the BIE)
PARAM.blockType = 1;                                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;                  %background flow
PARAM.deflationBlock = 1;

printToScreenStokes(PARAM);

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

figure
fx = yStokes(1:2:end-2);
fy = yStokes(2:2:end-1);
plot(Xsing,fx)
hold on
grid on
plot(Xsing,fy)
xlabel('x')
%legend('Analytical','Numerical','Location','Best')
ylabel('f_x,f_y')

U0 = yStokes(end);
U0an = 0.5/3;
errU = abs(U0an-U0)/abs(U0an);

display(['errU=' num2str(errU)])







