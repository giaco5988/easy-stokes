%compute velocity normal to the interface, the interface is reconstructed
%starting form modes coefficients

function [uMode,nx,ny,xGrid,yGrid,u] = NormalVelocityVolumeModes_BEM_curv(perturbMode,xBase,yBase,V0,xcm,PARAM)

    % numerics parameter
    PARAM.typeBCstokes = 2;
    PARAM.panels = [1];
    PARAM.cfunction = 0;
    PARAM.STstokes = 1;
    PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
    PARAM.panelType = [2];
    PARAM.blockType = [2];
    PARAM.addFlow = 2;
    PARAM.orderGeometryStokes = 1;
    PARAM.orderVariableStokes = 1;
    PARAM.deflationBlock = [0];
    PARAM.repulsiveForces = [0];
    PARAM.SPlinesType = [2];
    PARAM.stressBC = {1};

    % normal vector
    [nx,ny] = computeNormalVector(xBase',yBase',PARAM.orderVariableStokes(1),PARAM.orderGeometryStokes(1),PARAM.SPlinesType(1));

    %compute rho in symmetry axis
    if PARAM.dropFrame==0
        fVolume = @(firstMode) ModifyVolumeXYmodes(xBase,yBase,nx,ny,firstMode,V0,perturbMode);
    elseif PARAM.dropFrame==1
        fVolume = @(firstMode) ModifyVolumeXYmodesCM(xBase,yBase,nx,ny,firstMode,V0,perturbMode,xcm);
    end
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    firstMode = fsolve(fVolume,0,options);
    
    %build shape from mode
    perturbMode = [firstMode; perturbMode];
    perturb = chebcoeffs2chebvals(perturbMode);
    perturbCheb = chebfun(perturb,[0 1]);
    perturb = perturbCheb(linspace(0,1,numel(xBase)));
    
    %perturb respect to base shape
    xGrid = xBase' + nx.*perturb;
    yGrid = yBase' + ny.*perturb;
    
    %compute solution
    [y,~,~,nnx,nny] = BEM_Stokes({xGrid},{yGrid},PARAM);
    
    %normal velocity
    ux = y(1:2:end-1);  uy = y(2:2:end);
    nx = nnx{1}';   ny = nny{1}';
    u = nx.*ux + ny.*uy;
    
    % in drop frame
    if PARAM.dropFrame==1
        uDrop = DropVelocityAxis(xGrid,yGrid,u);
        u = nx.*(ux-uDrop) + ny.*uy;
    end
    
    %compute velocity modes coeff
    uCheb = chebfun(u,[0 1],'equi');
    uMode = chebcoeffs(uCheb);
    
    %take minus 1 velocty coefficient
    if numel(uMode)<PARAM.dealiasing
        uMode = [uMode(2:end); zeros(PARAM.dealiasing-numel(uMode),1)]; 
    else
        uMode = uMode(2:PARAM.dealiasing);
    end

end