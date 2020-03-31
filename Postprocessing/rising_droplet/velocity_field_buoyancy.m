%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p,X,Y,xxx,yyy] = velocity_field_buoyancy(a,b,solution,visc,capillary,xx1,xx2,yy,resx,resy,Replace)

    %set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

    %figure out center of mass velocity
    xcm = center_mass(a,b);

    %mesh finess per unit length
    radial = round(resy*yy);
    axial = round(resx*abs(xx1-xx2));

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
    U = (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %compute pressure
    %p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    p = -(PX*df_x'+PY*df_y')./InOut/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    
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

    %plot velocity field
%     figure
%     %subplot(1,3,1)
%     hold on
%     if drop_frame==1
%     quiver(Y,-X,V,-U-v_centermass,'b')
%     else
%     quiver(Y,-X,V,-U,'b')
%     end
%     axis equal
%     axis([-yy yy -xx2-0.05 -xx1])
%     %plot(Y0,-X0,'or')
%     %plot(b,-a,'gx')
%     if drop_frame==1
%     quiver(-Y,-X,-V,-U-v_centermass,'b')
%     else
%     quiver(-Y,-X,-V,-U,'b')
%     end
%     plot(yyy,-xxx,'k-',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
%     %quiver(b,-a,uy_sol',-ux_sol'-v_centermass,'r')
%     hold off
%     ylabel('axial direction')
%     xlabel('radial direction')
%     title(['flow field at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
%     
%     %plot pressure
% %     figure
% %     %subplot(1,3,2)
% %     contourf(Y,-X,p,200,'LineStyle','none')
% %     axis equal
% %     axis([-yy yy xx1 xx2])
% %     hold on
% %     contourf(-Y,-X,p,200,'LineStyle','none')
% %     plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
% %     hold off
% %     colorbar
% %     xlabel('axial direction')
% %     ylabel('radial direction')
% %     title(['Pressure for \lambda=' num2str(visc) ', Ca=' num2str(capillary)])
%     
%     %plot vorticity
%     figure
%     contourf(Y,-X,-omega,200,'LineStyle','none')
%     axis equal
%     axis([-yy yy xx1 xx2])
%     hold on
%     %plot(Y0,X0,'or')
%     contourf(-Y,-X,-omega,200,'LineStyle','none')
%     plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
%     hold off
%     colorbar
%     ylabel('axial direction')
%     xlabel('radial direction')
%     title(['Vorticity at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
%     
%     figure
%     plot(yyy,-xxx,'k',-yyy,-xxx,'k','LineWidth',2)
%     axis equal
%     axis([-1.5 1.5 -4.5 1.5])
%     title(['t=' num2str(ite*dt*checkpoint)])
    
end