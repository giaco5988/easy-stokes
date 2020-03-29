%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
L = 1;                    % motor lenght
h = 0.05;                    % wall thickness
manyTheta = linspace(pi/64,pi-pi/64,100);               % inclinitation of the cone

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
PARAM.deflationBlock = 1;

%geometry parameters
nPerLenght = 100;
PARAM.panels = 3;
PARAM.typePanel = [0 1 0];                % 0 is a straight line, 1 ia an arc
PARAM.typeBClaplace = [2 2 2];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0];     % 0 is straight, 1 is curved (spline)

%function profile for BCs
PARAM.fluxBC{1} = -1;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;

%filename
PARAM.filename = ['closedConeParametric_f1=' num2str(PARAM.fluxBC{1}) '_f2=' num2str(PARAM.fluxBC{2}) '_f3=' num2str(PARAM.fluxBC{3}) '_n=' num2str(nPerLenght) '_xi=' num2str(L) '_h=' num2str(h) '.mat'];

%initialize
U0 = zeros(numel(manyTheta),1);
Snum = zeros(numel(manyTheta),1);
solutionStokes = cell(numel(manyTheta),1);
solutionLaplace = cell(numel(manyTheta),1);
xxx = cell(numel(manyTheta),1);
yyy = cell(numel(manyTheta),1);

for i = 1:numel(manyTheta)
    
    theta = manyTheta(i);
    PARAM.n = [round((L-h/2/tan(theta))*nPerLenght) round(h*pi/2*nPerLenght) round((L+h/2/tan(theta))*nPerLenght)];                  % number of element per panel
    
    %print to screen
    printToScreenLaplace(PARAM)
    printToScreenStokes(PARAM)
    
    %build geometry
    xC = 0;   yC = 0;
    [x1,y1] = buildStraightLine(h/2/tan(theta),-h/2,L,-h/2,PARAM,1);
    [x2,y2] = buildArc(h/2,L,0,[-pi/2 pi/2],PARAM,2);
    [x3,y3] = buildStraightLine(L,h/2,-h/2/tan(theta),h/2,PARAM,3);
    x{1} = x1;  x{2} = x2;  x{3} = x3;
    y{1} = y1;  y{2} = y2;  y{3} = y3;

    %solid body rotation
    [x,y] = rigidBodyRotation(x,y,xC,yC,theta,1:3);

    %finite differences
    D1{1} = finiteDifference1D(PARAM.n(1),[2 0],1);
    D1{2} = finiteDifference1D(PARAM.n(2),[2 0],1);
    D1{3} = finiteDifference1D(PARAM.n(3),[2 0],1);

    %solve laplace equation
    yLaplace = BEM_Laplace(x,y,PARAM);

    %compute slip
    [vSlip1,l1] = computeSlipVelPhoretic(x{1},y{1},yLaplace(1:PARAM.n(1)),PARAM,1,D1{1});
    [vSlip2,l2] = computeSlipVelPhoretic(x{2},y{2},yLaplace(PARAM.n(1)+1:sum(PARAM.n(1:2))),PARAM,2,D1{2});
    [vSlip3,l3] = computeSlipVelPhoretic(x{3},y{3},yLaplace(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{3});

    %function profile for BCs
    PARAM.velBC{1} = -vSlip1*PARAM.fluxBC{1};
    PARAM.velBC{2} = -vSlip2*PARAM.fluxBC{2};
    PARAM.velBC{3} = -vSlip3*PARAM.fluxBC{3};

    %solve Stokes equation
    [yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);
    
    %save solutions
    solutionLaplace{i} = yLaplace;
    solutionStokes{i} = yStokes;
    xxx{i} = x;
    yyy{i} = y;

    %compute stresslet numerically
    Snum(i) = computeStressLetAxisNumerical(x,y,yStokes(1:end-1),PARAM);

    U0(i) = yStokes(end);
    display(['U0=' num2str(U0(i))])
    display(['S=' num2str(Snum(i))])

end

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

%plot velocity
figure
plot(manyTheta,U0)
grid on
xlabel('theta')
ylabel('U')
title('Motor velocity')

%plot stresslet
figure
plot(manyTheta,Snum)
grid on
xlabel('theta')
ylabel('S')
title('Stresslet intensity')









