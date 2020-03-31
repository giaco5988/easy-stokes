%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p,X,Y] = velocity_field_sphere(a,b,solution,xx1,xx2,yy,resx,resy,PARAM)

    set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

    visc = PARAM.visc;
    
    %mesh finess per unit length
    radial = round(resy*yy);
    axial = round(resx*abs(xx1-xx2));

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
    dfx_sol = solution(1:2:end-1);
    dfy_sol = solution(2:2:end);

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(0,yy,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = (InOut>pi)*0 + (InOut<pi);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_visu(a,b,X0,Y0);
    [PX,PY] = computeP_visu_const_StokesAX(a,b,X0,Y0);
    
    %compute the integral to obtain the velocities
    U = -PARAM.sphere_vel + (-GXX*dfx_sol - GXY*dfy_sol)/8/pi;
    V = (-GYX*dfx_sol - GYY*dfy_sol)/8/pi;
    
    %compute the integral to obtain the pressure
    p = -(PX*dfx_sol+PY*dfy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    %[X0,Y0,U,V,InOut] = close_interface_visu(a,b,X0,Y0,U,V,InOut,solution,1);
    
    %reshape in the grid shape
    X = reshape(X0,radial,axial);
    Y = reshape(Y0,radial,axial);
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    p = reshape(p,radial,axial);
    
    %moltiply time 1/lambda the points inside the droplet
    %InOut = reshape(InOut,radial,axial);
    %U = InOut.*U;
    %V = InOut.*V;
    
    %compute vorticity
    %omega = vorticity(axial,radial,U,V,X,Y);

    %plot velocity field
    figure
    quiver(X,Y,U,V)
    axis equal
    axis([-yy yy -xx2-0.05 -xx1])
    hold on
    %plot(-b,-a,'or')
    quiver(X,-Y,U,-V,'b')
    plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
    hold off
    ylabel('axial direction')
    xlabel('radial direction')
    title(['flow field at \lambda=' num2str(visc) ' and U=' num2str(PARAM.sphere_vel)])
    
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
    title(['Pressure for \lambda=' num2str(visc) ' and U=' num2str(PARAM.sphere_vel)])
    
    %plot vorticity
%     figure
%     contourf(Y,-X,-omega,200,'LineStyle','none')
%     axis equal
%     hold on
%     contourf(-Y,-X,-omega,200,'LineStyle','none')
%     plot(yyy,-xxx,'k',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
%     hold off
%     colorbar
%     ylabel('axial direction')
%     xlabel('radial direction')
%     title(['Vorticity at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(capillary) ', \delta=' num2str(delta)])
    
    
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