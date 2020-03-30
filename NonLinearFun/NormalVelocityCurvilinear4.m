%compute velocity normal to the interface

function [u,K1,K2,nx,ny] = NormalVelocityCurvilinear4(dn,xBase,yBase,lambda,PARAM)

    %compute first derivative in the curvilinear parameter
    [~,bx,cx,dx,~,by,cy,dy] = spline_symmetric(xBase',yBase');
    xpBase = derSplines(bx,cx,dx)';
    ypBase = derSplines(by,cy,dy)';

    %normal vector
    [nxBase,nyBase] = normalVectorSplineSymmetric(xpBase,ypBase);
    
    %by symmetry
    nxBase([1 end]) = [1 -1];
    nyBase([1 end]) = [0 0];

    %cartesian coordiantes
    x = xBase+dn.*nxBase;
    y = yBase+dn.*nyBase;
    
    %compute solution
%     [sol,nx,ny,K1,K2] = bemNewtonExtensUp(x,y,elem,lambda,capillary,PARAM);
    PARAM.typeBCstokes = 2;           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
    PARAM.orderVariableStokes = 1;    % 0 is constant on the elmennt, 1 is linear on the element
    PARAM.orderGeometryStokes = 1;    % 0 is straight, 1 is curved (spline)
    PARAM.SPlinesType = 2;            % 1 is natural splines, 2 is symmetric on the axis
    PARAM.panelType = 2;              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    PARAM.panelInflate = 0;           % 0 is not inflating, 1 is inflating
    PARAM.blockType = 2;                    % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
    PARAM.deflationBlock = 0;
    PARAM.deflationConstant = 0;
    PARAM.addFlow = 0;
    PARAM.Qsource = 0;
    PARAM.cfunction = 0;
    PARAM.STstokes = 1;
    PARAM.kernelFreeSpace = 1;
    PARAM.posWall = nan;
    PARAM.visc = lambda;
    PARAM.stressBC{1} = 1;
    PARAM.panels = [1];
    PARAM.repulsiveForces = [0];
    PARAM.addFlow = 2;
    [sol,~,~,nx,ny] = BEM_Stokes({x'},{y'},PARAM);
    nx = nx{1}'; ny = ny{1}';
    [K1,K2] = computeCurvatureSplines(x',y',1);
    
    %axial and radial
    ux = sol(1:2:end-1);  uy = sol(2:2:end);
    
    %compute normal velocity
    if PARAM.dropFrame==0
        u = nx.*ux + ny.*uy;
    elseif PARAM.dropFrame==1
        u = NormalVelocityDropFrame(x',y',ux,uy);
    end

end