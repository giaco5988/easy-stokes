%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,X,Y] = velocity_field_uniform(a,b,solution,visc,capillary,xx1,xx2,yy,resx,resy,dt,checkpoint,delta,PARAM)

    set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

    %figure out center of mass velocity
    %v_centermass = -(xcm-xcm_old)/dt/checkpoint;    %minus beacuse I want a boyant droplet

    %visc2 = 1;
    %visc1 = visc*visc2;

    %mesh finess per unit length
    radial = round(resy*yy);
    axial = round(resx*abs(xx1-xx2));

    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
    
    %compute the versor normal tp the node (N) and to the element (r)
    %[N,r] =normal_versor_clean(a,b,q-1);
    N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
  
   %%%%%%%%%%%%%%%%%%% COMPUTE CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%%%%%%%

   K1 = curv_spline2(bx,by,cx,cy,dx,dy);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   %those because the points are on the axis
   N(:,1) = [1; 0];
   N(:,end) = [-1; 0];
  
   %second component of the curvature
   K2 = N(2,:)./b;
   K2(1) = 0;
   K2(end) = 0;
  
   %this because on the axis the I cannot use the definition of before
   K1(1) = 2*K1(1);
   K1(end) = 2*K1(end);
  
   %sum the two component of the curvature
   K = K1+K2;

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
    cap_forces = K/capillary;
    df = cap_forces;
    [df_x,df_y] = stress_diff(df,N,numel(a));
    
    %BC velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);

    %create the grid in cartesina coordinate
    x = linspace(xx1,xx2,axial);
    y = linspace(0,yy,radial);
    
    %create the grid in cylindrical coordinate
    %r = linspace(0,xx2,axial);
    %theta = linspace(0,2*pi,radial);

    [X,Y] = meshgrid(x,y);
    %[R,THETA] = meshgrid(r,theta);
    %X = R.*cos(THETA);
    %Y = R.*sin(THETA);
    
    X0 = X(:);
    Y0 = Y(:);

    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = (InOut>pi)/visc + (InOut<pi);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    [PX,PY,PI1,PI2] = computeP_AX_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    
    %compute the integral to obtain the velocities
    U = PARAM.U + (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %compute pressure
    p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    [X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,1);
    
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
    %subplot(1,3,1)
    hold on
    quiver(Y,X,-V,U-PARAM.U,'b')
    axis equal
    axis([-yy yy -xx2-0.05 -xx1])
    %plot(Y0,-X0,'or')
    quiver(-Y,X,V,U-PARAM.U,'b')
    plot(yyy,-xxx,'k-',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
    %quiver(b,-a,uy_sol',-ux_sol'-v_centermass,'r')
    hold off
    ylabel('axial direction')
    xlabel('radial direction')
    title(['flow field for \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
    
    %plot pressure
    figure
    %subplot(1,3,2)
    contourf(Y,-X,p,200,'LineStyle','none')
    axis equal
    axis([-yy yy xx1 xx2])
    hold on
    contourf(-Y,-X,p,200,'LineStyle','none')
    plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    xlabel('axial direction')
    ylabel('radial direction')
    title(['Pressure for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
    
    %plot vorticity
    figure
    %subplot(1,3,3)
    contourf(Y,-X,-omega,200,'LineStyle','none')
    axis equal
    axis([-yy yy xx1 xx2])
    hold on
    contourf(-Y,-X,-omega,200,'LineStyle','none')
    plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
    hold off
    colorbar
    ylabel('axial direction')
    xlabel('radial direction')
    title(['Vorticity at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
    
    
%     figure
%     contourf(Y,-X,sqrt((U+v_centermass).^2+V.*V),200,'LineStyle','none')
%     colorbar
%     axis equal
%     hold on
%     contourf(-Y,-X,sqrt((U+v_centermass).^2+V.*V),200,'LineStyle','none')
%     plot(yyy,-xxx,'r',-yyy,-xxx,'r','MarkerSize',2,'LineWidth',2)
%     hold off
%     ylabel('axial direction')
%     xlabel('radial direction')
%     title(['Velocity magnitude at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
    
end