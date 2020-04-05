%compute velocity normal to the interface, the interface is reconstructed
%starting form modes coefficients

function [uMode,nx,ny,xGrid,yGrid,u,xBase,yBase] = NormalVelocityPipeModesBEM(perturbMode,xBase,yBase,V0,xcm,PARAM)

    %% fixed parameters
    dropFrame = 1;
    PARAM.cfunction = 0;
    PARAM.visc = [0 0 0 PARAM.visc];
    PARAM.repulsiveForces = [0 0];
    
    %% numerics parameter
    L = 20;
    H = 1;
    nPerLenght = 10;
    PARAM.panels = [3 1];                   % number of panels per block (one block is a closed line, 1 panel is each boundary where BCs are applied)
    PARAM.rotate = [0 0 0 0]/180*pi;        % geometry element rotation        
    PARAM.n = [round(H*nPerLenght)
        round(L*nPerLenght)
        round(H*nPerLenght)
        numel(xBase)-1];                        % number of element per panel
    PARAM.orderVariableStokes = [0 0 0 1];    % 0 is constant on the elmennt, 1 is linear on the element
    PARAM.orderGeometryStokes = [0 0 0 1];    % 0 is straight, 1 is curved (spline)
    PARAM.SPlinesType = [0 0 0 2];            % 1 is natural splines, 2 is symmetric on the axis
    PARAM.panelType = [0 0 0 2];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    PARAM.panelInflate = [0 0 0 0];           % 0 is not inflating, 1 is inflating
    PARAM.blockType = [0 2];                  % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    PARAM.deflationBlock = [0 0];             % whether to apply deflation which is necessary for 0 viscoisty ratio and moving rigid objects (no by default)
    PARAM.deflationConstant = [0 4*pi];       % constants associated to deflation
    PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
    PARAM.STstokes = 1;
    
    %% Geometry parameters
    PARAM.geometryPanel = [0 0 0 2];        % 0 is a straight line, 1 ia an arc, 2 is a defromable object
    PARAM.xStart = [-L/2 -L/2 L/2 nan];     % x starting point for the straight lines
    PARAM.xEnd = [-L/2 L/2 L/2 nan];        % x ending point for the straight lines
    PARAM.yStart = [0 H H nan];             % y starting point for the straight lines
    PARAM.yEnd = [H H 0 nan];               % y ending point for the straight lines
    PARAM.thetaStart = [nan nan nan 0];     % theta starting point for the arc
    PARAM.thetaEnd = [nan nan nan pi];      % theta starting point for the arc
    PARAM.rArc = [nan nan nan PARAM.alpha]; % theta starting point for the arc
    PARAM.x0_Circle = [nan nan nan xcm];    % center of the circle (axial direction)
    PARAM.y0_Circle = [nan nan nan 0];      % center of the circle (radial direction)
    PARAM.xCrotate = [0 0 0 0];             % rotation for circles
    PARAM.yCrotate = [0 0 0 0];             % rotation for circles
    PARAM.ellipseShape = [0 0 0 1];         % rotation for circles
    PARAM.D = [0 0 0 1];
    [xGeometry,yGeometry,PARAM] = buildGeometryPanelsGeneral(PARAM);
    
    %% BCs
    PARAM.typeBCstokes = [8 0 0 2];
    PARAM.addFlow = 0;                        % add background flow (for example extensional flow)
    PARAM.velBCaxial = {nan 0 2*PARAM.Ca*(1-(yGeometry{3}(1:end-1)+yGeometry{3}(2:end))/2.^2) nan};
    PARAM.velBCradial = {0 0 0 nan};
    PARAM.stressBC = {0 nan nan 1};

    % compute normal vector
    [nx,ny] = computeNormalVector(xBase',yBase',PARAM.orderVariableStokes(4),PARAM.orderGeometryStokes(4),PARAM.SPlinesType(4));

    %compute two complementary modes, fue to volume and center of mass
    %constraint
    fVolume = @(firstMode) ModifyVolumeXYmodesCM(xBase,yBase,nx,ny,firstMode,V0,perturbMode,xcm);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    firstMode = fsolve(fVolume,[0; 0],options);
    
    % build shape from mode
    perturbMode = [firstMode; perturbMode];
    perturb = chebcoeffs2chebvals(perturbMode);
    perturbCheb = chebfun(perturb,[0 1]);
    perturb = perturbCheb(linspace(0,1,numel(xBase)));
    
    % perturb respect to base shape
    xGrid = xBase' + nx.*perturb;
    yGrid = yBase' + ny.*perturb;
    
    % compute solution
    xGeometry{4} = xGrid;   yGeometry{4} = yGrid;
    [y,~,~,nnx,nny] = BEM_Stokes(xGeometry,yGeometry,PARAM);
    
    % normal velocity
    ux = y(2*sum(PARAM.n(1:3))+1:2:end-1);  uy = y(2*sum(PARAM.n(1:3))+2:2:end);
    nx = nnx{4}';   ny = nny{4}';
    u = nx.*ux + ny.*uy;
    
    % in drop frame
    if dropFrame==1
        uNormal = @(uDrop) nx.*(ux-uDrop) + ny.*uy;
        fVelDrop = @(Udrop) DropVelocityAxisNewton(xGrid,yGrid,uNormal,Udrop);
        options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
        VelDrop = fsolve(fVelDrop,DropVelocityAxis(xGrid,yGrid,u),options);
        u = uNormal(VelDrop);
    end
    
    % compute velocity modes coeff
    uCheb = chebfun(u,[0 1],'equi');
    uMode = chebcoeffs(uCheb);
    
    % run checks
    if numel(firstMode) ~= 2
            error('Should have two compensating modes when in drop frame')
    
    end
    
    % take minus 2 velocity coefficient
    if numel(uMode)<PARAM.dealiasing
        error('not implemented')
    else
        uMode = uMode(3:PARAM.dealiasing+2);
    end
    
    if numel(uMode) ~= numel(perturbMode)-numel(firstMode)
        error('Number of fval has to be the same as DOF')
    end

end