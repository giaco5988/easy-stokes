%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,X,Y] = velocity_field_2drop_poiseuille(a,b,solution,xx1,xx2,yy,resx,resy,ite,PARAM)

    set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
    
    visc = PARAM.visc;
    capillary1 = PARAM.Ca1;
    capillary2 = PARAM.Ca2;
    dt = PARAM.deltaT;
    checkpoint = PARAM.checkpoint;
    
    %compute gamma with the CURRENT DEFINITION OF CA
    U = PARAM.Q/pi/PARAM.R^2;
    gamma1 = U*PARAM.visc2/capillary1;
    gamma2 = U*PARAM.visc2/capillary2;

    solution = solution(:,ite);
    a_before = a(:,ite-2)';
    b_before = b(:,ite-2)';
    a = a(:,ite-1)';
    b = b(:,ite-1)';
    
    %figure out center of mass velocity
    xcm_old = center_mass(a_before,b_before);
    xcm = center_mass(a,b);
    v_centermass = (xcm-xcm_old)/dt/checkpoint;
    
    %center of mass
    a = a-xcm;
    
    q = PARAM.q;
    %p = PARAM.p;
    
    %interface 1 and 2
    a1 = a(1:q+1);
    b1 = b(1:q+1);
    a2 = a(q+2:end);
    b2 = b(q+2:end);

    %mesh finess per unit length
    radial = round(resy*yy);
    axial = round(resx*abs(xx1-xx2));

    %compute the spline coeff
    [ax1, bx1, cx1, dx1, ay1, by1, cy1, dy1] = spline_symmetric (a1, b1);
    [ax2, bx2, cx2, dx2, ay2, by2, cy2, dy2] = spline_symmetric (a2, b2);

    %compute splines coordinates
    t = 0:0.1:0.9;
    
    %INTERFACE 1
    ttt = repmat(t,1,numel(ax1));
    axxx1 = reshape(repmat(ax1,numel(t),1),1,numel(ax1)*numel(t));
    bxxx1 = reshape(repmat(bx1,numel(t),1),1,numel(bx1)*numel(t));
    cxxx1 = reshape(repmat(cx1,numel(t),1),1,numel(cx1)*numel(t));
    dxxx1 = reshape(repmat(dx1,numel(t),1),1,numel(dx1)*numel(t));
    ayyy1 = reshape(repmat(ay1,numel(t),1),1,numel(ay1)*numel(t));
    byyy1 = reshape(repmat(by1,numel(t),1),1,numel(by1)*numel(t));
    cyyy1 = reshape(repmat(cy1,numel(t),1),1,numel(cy1)*numel(t));
    dyyy1 = reshape(repmat(dy1,numel(t),1),1,numel(dy1)*numel(t));
  
    %splines coordinates
    xxx1 = [axxx1+bxxx1.*ttt+cxxx1.*ttt.^2+dxxx1.*ttt.^3 a1(end)];
    yyy1 = [ayyy1+byyy1.*ttt+cyyy1.*ttt.^2+dyyy1.*ttt.^3 b1(end)];
    
    %INTERFACE 2
    ttt = repmat(t,1,numel(ax2));
    axxx2 = reshape(repmat(ax2,numel(t),1),1,numel(ax2)*numel(t));
    bxxx2 = reshape(repmat(bx2,numel(t),1),1,numel(bx2)*numel(t));
    cxxx2 = reshape(repmat(cx2,numel(t),1),1,numel(cx2)*numel(t));
    dxxx2 = reshape(repmat(dx2,numel(t),1),1,numel(dx2)*numel(t));
    ayyy2 = reshape(repmat(ay2,numel(t),1),1,numel(ay2)*numel(t));
    byyy2 = reshape(repmat(by2,numel(t),1),1,numel(by2)*numel(t));
    cyyy2 = reshape(repmat(cy2,numel(t),1),1,numel(cy2)*numel(t));
    dyyy2 = reshape(repmat(dy2,numel(t),1),1,numel(dy2)*numel(t));
  
    %splines coordinates
    xxx2 = [axxx2+bxxx2.*ttt+cxxx2.*ttt.^2+dxxx2.*ttt.^3 a2(end)];
    yyy2 = [ayyy2+byyy2.*ttt+cyyy2.*ttt.^2+dyyy2.*ttt.^3 b2(end)];
  
    %BC stresses
    [df_x1,df_y1] = df_surf_tens_spline_symm(a1,b1,gamma1);
    [df_x2,df_y2] = df_surf_tens_spline_symm(a2,b2,gamma2);
    
    %merge stresses difference
    df_x = [df_x1 df_x2];
    df_y = [df_y1 df_y2];
    
    %BC velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(0,yy,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %figure out if I'm inside or outside the droplet2
    InOut1 = FigInOut(a1,b1,X0,Y0);
    InOut2 = FigInOut(a2,b2,X0,Y0);
    
    InOut = InOut1 + InOut2;
    InOut = (InOut>pi)/visc + (InOut<pi);

    %compute the necessary Green's function INTEGRATION ON INTERFACE 1
    [GXX1,GXY1,GYX1,GYY1,A111,A121,A211,A221] = computeGT_spline_visu(ax1,ay1,bx1,by1,cx1,cy1,dx1,dy1,X0',Y0');
    
    %compute the necessary Green's function INTEGRATION ON INTERFACE 2
    [GXX2,GXY2,GYX2,GYY2,A112,A122,A212,A222] = computeGT_spline_visu(ax2,ay2,bx2,by2,cx2,cy2,dx2,dy2,X0',Y0');
    
    %concatenate
    GXX = [GXX1 GXX2];
    GXY = [GXY1 GXY2];
    GYX = [GYX1 GYX2];
    GYY = [GYY1 GYY2];
    A11 = [A111 A112];
    A12 = [A121 A122];
    A21 = [A211 A212];
    A22 = [A221 A222];
    
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
    
    [PX1,PY1,PI11,PI21] = computeP_AX_spline_visu(ax1,ay1,bx1,by1,cx1,cy1,dx1,dy1,X0',Y0');
    [PX2,PY2,PI12,PI22] = computeP_AX_spline_visu(ax2,ay2,bx2,by2,cx2,cy2,dx2,dy2,X0',Y0');
    %[PX1,PY1] = computeP_visu_const_StokesAX(a1,b1,X0,Y0);
    %[PX2,PY2] = computeP_visu_const_StokesAX(a2,b2,X0,Y0);
    
    PX = [PX1 PX2];
    PY = [PY1 PY2];
    PI1 = [PI11 PI12];
    PI2 = [PI21 PI22];
    
    %compute the integral to obtain the velocities with undelying
    %velocities
    u_x = poiseuille_flow(Y0,PARAM.Q,PARAM.R);
    U = u_x + (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %compute the integral to obtain the pressure
    %df_x = (df_x(1:end-1)+df_x(2:end))/2;
    %df_y = (df_y(1:end-1)+df_y(2:end))/2;
    %p = -(PX*df_x'+PY*df_y')/8/pi;
    p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    %[X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,1);
    
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
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    quiver(X,-Y,U-v_centermass,-V,'b')
    plot(xxx1,yyy1,'k',xxx1,-yyy1,'k','MarkerSize',2,'LineWidth',2)
    plot(xxx2,yyy2,'k',xxx2,-yyy2,'k','MarkerSize',2,'LineWidth',2)
    hold off
    xlabel('axial direction')
    ylabel('radial direction')
    title(['flow field for \lambda=' num2str(visc) ', Ca1=' num2str(capillary1) ' Ca2=' num2str(capillary2)])
    
    %plot pressure
    figure
    contourf(X,Y,p,200,'LineStyle','none')
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    contourf(X,-Y,p,200,'LineStyle','none')
    plot(xxx1,yyy1,'k',xxx1,-yyy1,'k','MarkerSize',2,'LineWidth',2)
    plot(xxx2,yyy2,'k',xxx2,-yyy2,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    xlabel('axial direction')
    ylabel('radial direction')
    title(['Pressure for \lambda=' num2str(visc) ', Ca1=' num2str(capillary1) ' Ca2=' num2str(capillary2)])
    
    %plot vorticity
    figure
    contourf(X,Y,omega,200,'LineStyle','none')
    axis equal
    axis([xx1 xx2 -yy yy])
    hold on
    contourf(X,-Y,omega,200,'LineStyle','none')
    quiver(X,-Y,U-v_centermass,-V,'k')
    quiver(X,Y,U-v_centermass,V,'k')
    plot(xxx1,yyy1,'k',xxx1,-yyy1,'k','MarkerSize',2,'LineWidth',2)
    plot(xxx2,yyy2,'k',xxx2,-yyy2,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    xlabel('axial direction')
    ylabel('radial direction')
    title(['Vorticity for \lambda=' num2str(visc) ', Ca1=' num2str(capillary1) ' Ca2=' num2str(capillary2)])

end