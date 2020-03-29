%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
%close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%parameter
manyTheta = linspace(0,pi,101);
%manyTheta = linspace(0,pi,6);
U0 = zeros(numel(manyTheta),1);

for i = 1:numel(manyTheta)
    
    display([num2str(i) ' of ' num2str(numel(manyTheta)) ' theta=' num2str(manyTheta(i))])

    %physical parameters
    r = 1;                      % average cone radius
    L = 1.8;                     % motor lenght
    h = 0.2;                    % wall thickness
    theta = manyTheta(i);       % inclinitation of the cone
    PARAM.phi = [0 0];          % concentration
    PARAM.flux = [0 0 1 0];     % BC for flux

    %options
    PARAM.cfunction = 0;
    PARAM.STlaplace = 1;
    PARAM.STstokes = 1;
    plotField = 0;

    %geometry parameters
    nPerLenght = 50;
    PARAM.panels = 4;
    PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel
    PARAM.typePanel = [0 1 0 1];                % 0 is a straight line, 1 ia an arc
    PARAM.typeBClaplace = [2 2 2 2];            % 1 is prescribed concentration, 2 is prescibed flux
    PARAM.orderVariableLaplace = [0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
    PARAM.orderGeometryLaplace = [0 0 0 0];     % 0 is straight, 1 is curved (spline)

    %function profile for BCs
    PARAM.fluxBC{1} = 1;
    PARAM.fluxBC{2} = 1;
    PARAM.fluxBC{3} = 1;
    PARAM.fluxBC{4} = 1;
    
    %numerics parameters for Stokes
    PARAM.typeBCstokes = [3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
    PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
    PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
    PARAM.panelType = [1 1 1 1];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    
    if i==1
        %print to screen
        printToScreenLaplace(PARAM)
        printToScreenStokes(PARAM)
    end

    %build geometry
    xC = L/2;   yC = r;
    [x1,y1] = buildStraightLine(L,r+h/2,0,r+h/2,PARAM,1);
    [x2,y2] = buildArc(h/2,0,r,[pi/2 3*pi/2],PARAM,2);
    [x3,y3] = buildStraightLine(0,r-h/2,L,r-h/2,PARAM,3);
    [x4,y4] = buildArc(h/2,L,r,[-pi/2 pi/2],PARAM,4);
    x{1} = x1;  x{2} = x2;  x{3} = x3;  x{4} = x4;
    y{1} = y1;  y{2} = y2;  y{3} = y3;  y{4} = y4;

    %solid body rotation
    [x,y] = rigidBodyRotation(x,y,xC,yC,theta,4);

    %finite differences
    D1{1} = finiteDifference1D(PARAM.n(1),[2 0],1);
    D1{2} = finiteDifference1D(PARAM.n(2),[2 0],1);
    D1{3} = finiteDifference1D(PARAM.n(3),[2 0],1);
    D1{4} = finiteDifference1D(PARAM.n(4),[2 0],1);
    
    %filename
    PARAM.filename = ['conePhoreticAllMob_f1=' num2str(PARAM.flux(1)) '_f2=' num2str(PARAM.flux(2)) '_f3=' num2str(PARAM.flux(3)) '_f4=' num2str(PARAM.flux(4)) '_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

    %solve laplace equation
    conc = BEM_Laplace(x,y,PARAM);

    %compute slip
    [vSlip1,l1] = computeSlipVelPhoretic(x{1},y{1},conc(+1:PARAM.n(1)),PARAM,1,D1{1});
    [vSlip2,l2] = computeSlipVelPhoretic(x{2},y{2},conc(sum(PARAM.n(1))+1:sum(PARAM.n(1:2))),PARAM,2,D1{2});
    [vSlip3,l3] = computeSlipVelPhoretic(x{3},y{3},conc(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{3});
    [vSlip4,l4] = computeSlipVelPhoretic(x{4},y{4},conc(sum(PARAM.n(1:3))+1:sum(PARAM.n(1:4))),PARAM,4,D1{4});

    %function profile for BCs
    PARAM.velBC{1} = vSlip1;
    PARAM.velBC{2} = vSlip2;
    PARAM.velBC{3} = vSlip3;
    PARAM.velBC{4} = vSlip4;
    
    %solve Stokes equation
    [yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

    %motor velocity
    U0(i) = yStokes(end);
    
    %save results
    cd(PARAM.res)
    save(PARAM.filename)
    cd(PARAM.here)

end
















