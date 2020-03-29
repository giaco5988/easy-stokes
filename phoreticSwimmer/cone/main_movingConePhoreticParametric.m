%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%parameter
manyTheta = linspace(1e-3,pi-1e-3,100);
L = 1;
nPerLenght = 100;            % elemnts per unit lenght
manyR = 0.25;                      % radial coordinate of pivoting point
h = 0.05;

%initialize
U0 = zeros(numel(manyTheta),numel(manyR));
Snum = zeros(numel(manyTheta),numel(manyR));
solutionStokes = cell(numel(manyTheta),numel(manyR));
solutionLaplace = cell(numel(manyTheta),numel(manyR));
xxx = cell(numel(manyTheta),numel(manyR));
yyy = cell(numel(manyTheta),numel(manyR));

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;

%geometry parameters
PARAM.panels = 4;                           % number of element per panel
PARAM.typePanel = [0 1 0 1];                % 0 is a straight line, 1 ia an arc
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel

%numerics parameters for Laplace
PARAM.typeBClaplace = [2 2 2 2];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0];     % 0 is straight, 1 is curved (spline)

%function profile for BCs Laplace
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
PARAM.deflationBlock = 1;

%filename
if numel(manyR)>1
    PARAM.filename = ['conePhoreticParametric_n=' num2str(nPerLenght) '_L=' num2str(L) '_h=' num2str(h) '_f1=' num2str(PARAM.fluxBC{1}) '_f2=' num2str(PARAM.fluxBC{2}) '_f3=' num2str(PARAM.fluxBC{3}) '_f4=' num2str(PARAM.fluxBC{4}) '.mat'];
elseif numel(manyR)==1
    PARAM.filename = ['conePhoreticParametric_n=' num2str(nPerLenght) '_L=' num2str(L) '_r=' num2str(manyR) '_h=' num2str(h) '_f1=' num2str(PARAM.fluxBC{1}) '_f2=' num2str(PARAM.fluxBC{2}) '_f3=' num2str(PARAM.fluxBC{3}) '_f4=' num2str(PARAM.fluxBC{4}) '.mat'];
end

%simulation settings
printToScreenLaplace(PARAM);
printToScreenStokes(PARAM);
    
for k = 1:numel(manyR)

    r = manyR(k);                     % height of pivoting point
    
    display([num2str(k) ' of ' num2str(numel(manyR)) ' r=' num2str(manyR(k))])
    
    for i = 1:numel(manyTheta)

        %physical parameters
        theta = manyTheta(i);       % inclinitation of the cone

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

        %solve laplace equation
        [yLaplace,~,~,nnn] = BEM_Laplace(x,y,PARAM);

        %decide velocity BC for Stokes
        for l = 1:4
            
            if PARAM.fluxBC{l}~=0
                %finite differences
                D1 = finiteDifference1D(PARAM.n(l),[2 0],1);
                
                if l==1
                    startMatrix = 1;
                else
                    startMatrix = 1+sum(nnn(1:l-1));
                end
                endMatrix = sum(nnn(1:l));
                
                %compute slip velocity
                vSlip = computeSlipVelPhoretic(x{l},y{l},yLaplace(startMatrix:endMatrix),PARAM,l,D1);

                %function profile for BCs
                PARAM.velBC{l} = vSlip;
                
            else
                
                PARAM.velBC{l} = 0;
        
            end
        
        end

        %solve Stokes equation
        [yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

        %motor velocity
        U0(i,k) = yStokes(end);
        
        %save solutions
        solutionLaplace{i,k} = yLaplace;
        solutionStokes{i,k} = yStokes;
        xxx{i,k} = x;
        yyy{i,k} = y;
        
        %stresslet
        Snum(i,k) = computeStressLetAxisNumerical(x,y,yStokes(1:end-1),PARAM);   

    end
end

%save results
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

if numel(manyR)>1
    
    %motor velocity
    figure
    contourf(manyTheta,manyR/L,U0')
    xlabel('\theta')
    ylabel('L/r')
    title('Motor Velocity')
    grid on
    colorbar

    %stresslet
    figure
    contourf(manyTheta,manyR/L,Snum')
    xlabel('\theta')
    ylabel('L/r')
    title('Stresslet intensity')
    grid on
    colorbar
    
elseif numel(manyR)==1
    
    %plot velocity
    figure
    plot(manyTheta,U0)
    grid on
    xlabel('\theta')
    ylabel('U')
    title('Motor velocity')

    %plot stresslet
    figure
    plot(manyTheta,Snum)
    grid on
    xlabel('\theta')
    ylabel('S')
    title('Stresslet intensity')

end















