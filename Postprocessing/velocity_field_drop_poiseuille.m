%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,X,Y] = velocity_field_drop_poiseuille(a,b,solution,xx1,xx2,yy,resx,resy,ite,PARAM)

    set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
    
    visc = PARAM.visc;
    capillary = PARAM.Ca;
    dt = PARAM.deltaT;
    checkpoint = PARAM.checkpoint;
    
    %compute gamma with the CURRENT DEFINITION OF CA
    U = PARAM.Q/pi/PARAM.R^2;
    gamma = U*PARAM.visc2/capillary;

    solution = solution(:,ite);
    a_before = a(:,ite-2)';
    b_before = b(:,ite-2)';
    a = a(:,ite-1)';
    b = b(:,ite-1)';

    %figure out center of mass velocity
    xcm_old = center_mass(a_before,b_before);
    xcm = center_mass(a,b);
    v_centermass = (xcm-xcm_old)/dt/checkpoint;

    %mesh finess per unit length
    radial = round(resy*yy);
    axial = round(resx*abs(xx1-xx2));

    %center of mass
    %xcm = center_mass(a,b);
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
    [df_x,df_y] = df_surf_tens_spline_symm(a,b,gamma);
    
    %BC velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(0,yy,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);
    
%     figure
%     plot(X0,Y0,'o')

    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = (InOut>pi)/visc + (InOut<pi);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    
    %compute the necessary Green's function with C functions
%     [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] =...
%             Stokes2DAxisSPlinesLinear(ax,bx,cx,dx,ay,by,cy,dy,X0',Y0');
%         
%     q = PARAM.q;
%     sing = numel(X0);
%         
%     GXX = reshape(GXX,sing,q+1);
%     GXY = reshape(GXY,sing,q+1);
%     GYX = reshape(GYX,sing,q+1);
%     GYY = reshape(GYY,sing,q+1);
%     TXXX = reshape(TXXX,sing,q+1);
%     TXXY = reshape(TXXY,sing,q+1);
%     TXYX = reshape(TXYX,sing,q+1);
%     TXYY = reshape(TXYY,sing,q+1);
%     TYXX = reshape(TYXX,sing,q+1);
%     TYXY = reshape(TYXY,sing,q+1);
%     TYYX = reshape(TYYX,sing,q+1);
%     TYYY = reshape(TYYY,sing,q+1);
%     
%     A11 = TXXX + TXXY;
%     A12 = TXYX + TXYY;
%     A21 = TYXX + TYXY;
%     A22 = TYYX + TYYY;
    
    [PX,PY,PI1,PI2] = computeP_AX_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    %[PX,PY] = computeP_visu_const_StokesAX(a,b,X0,Y0);
    
    %compute the integral to obtain the velocities with undelying
    %velocities
    u_x = poiseuille_flow(Y0,PARAM.Q,PARAM.R);
    U = u_x + (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %compute the integral to obtain the pressure
    p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    [X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,2);
    
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
    omega = vorticity(axial,radial,U,V,X,Y);

    %plot velocity field
    figure
    quiver(X,Y,U-v_centermass,V)
    %quiver(X,Y,U,V)
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    quiver(X,-Y,U-v_centermass,-V,'b')
    %quiver(a,b,ux_sol'-v_centermass,uy_sol','r')
    %quiver(a,-b,ux_sol'-v_centermass,-uy_sol','r')
    plot(xxx,yyy,'k',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
    hold off
    xlabel('axial direction')
    ylabel('radial direction')
    title(['flow field for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
    
    %plot pressure
    figure
    contourf(X,Y,p,200,'LineStyle','none')
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    contourf(X,-Y,p,200,'LineStyle','none')
    plot(xxx,yyy,'k',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    xlabel('axial direction')
    ylabel('radial direction')
    title(['Pressure for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
    
    %plot vorticity
    figure
    contourf(X,Y,omega,200,'LineStyle','none')
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    contourf(X,-Y,omega,200,'LineStyle','none')
    plot(xxx,yyy,'k',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    xlabel('axial direction')
    ylabel('radial direction')
    title(['Vorticity for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
    
%     figure
%     contourf(X,Y,sqrt((U-v_centermass).^2+V.*V),200,'LineStyle','none')
%     colorbar
%     axis equal
%     axis([xx1 xx2 -yy yy])
%     hold on
%     contourf(X,-Y,sqrt((U-v_centermass).^2+V.*V),200,'LineStyle','none')
%     plot(yyy,-xxx,'r',-yyy,-xxx,'r','MarkerSize',2,'LineWidth',2)
%     hold off
%     ylabel('axial direction')
%     xlabel('radial direction')
%     title(['Velocity magnitude for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
    
end