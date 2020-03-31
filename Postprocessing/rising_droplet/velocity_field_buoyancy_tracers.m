%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p,X,Y,xxx,yyy] = velocity_field_buoyancy_tracers(a,b,solution,visc,capillary,X,Y,Replace,dropFrame)

    %set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

    %figure out center of mass velocity
    xcm = center_mass(a,b);

    %center of mass
    a = a-xcm;

    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);

    %compute splines coordinates
    t = 0:0.1:0.9;
    ttt = repmat(t,1,numel(ax));
    axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
    bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
    cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
    dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
    ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
    byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
    cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
    dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
  
    %splines coordinates
    xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a(end)];
    yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b(end)];
  
    %BC stresses
    [df_x,df_y] = df_buoyancy_spline_symm(a,b,capillary,visc,0);
    
    %BC velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);
    
    X0 = X(:);
    Y0 = Y(:);
    
    temp = size(X);
    radial = temp(1); axial = temp(2);

    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = (InOut>pi)/visc + (InOut<pi);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    [PX,PY,PI1,PI2] = computeP_AX_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    
    %compute the integral to obtain the velocities
    U = (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %compute pressure
    p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    [X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,Replace);
    
    %reshape in the grid shape
    X = reshape(X0,radial,axial);
    Y = reshape(Y0,radial,axial);
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    p = reshape(p,radial,axial);
    
    %moltiply time 1/lambda the points inside the droplet
    InOut = reshape(InOut,radial,axial);
    U = InOut.*U;
    V = InOut.*V;
    
    %compute vorticity
    %omega = vorticity(axial,radial,U,V,X,Y);

    if dropFrame==1
        
        %compute and subtract droplet velocity
        Un = DropNormalVelocity(a,b,ux_sol,uy_sol);
        Vdrop = DropVelocityAxis(a,b,Un);
        U = U-Vdrop;
        
    end
    
end