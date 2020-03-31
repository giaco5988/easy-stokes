%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p,X,Y,xxx,yyy] = velocity_field_extensional(a,b,solution,visc,capillary,xx1,xx2,yy,resx,resy,Replace)

  %set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',2,'DefaultAxesFontName','Times New Roman');

  %visc2 = 1;
  %visc1 = visc*visc2;
  
  %mesh finess
  %radial = 5;
  %axial = 5;
  
  %mesh finess per unit length
  radial = round(resy*yy);
  axial = round(resx*abs(xx1-xx2));
    
  xcm = center_mass(a,b);
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
  
  %compute the versor normal tp the node (N) and to the element (r)
  %[N,r] =normal_versor_clean(a,b,q-1);
  N = normal_splines(bx,cx,dx,by,cy,dy);
  
  %those because the points are on the axis
  N(:,1) = [1; 0];
  N(:,end) = [-1; 0];
  
  %%%%%%%%%%%%%% COMPUTE MERIDIONAL CURVATURE WITH SPLINES %%%%%%%%%%%%%%%%

  K1 = curv_spline2(bx,by,cx,cy,dx,dy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%% COMPUTE AZIMUTHAL CURVATURE %%%%%%%%%%%%%%%%%%%%%%%%
  
  %second component of the curvature
  K2 = N(2,:)./b;
  
  %this because on the axis the I cannot use the definition of before
  K2(1) = K1(1);
  K2(end) = K1(end);
  
  %sum the two component of the curvature
  K = K1+K2;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE DELTA_F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  cap_forces = K;
  
  df = cap_forces;
  [df_x,df_y] = stress_diff(df,N,numel(a));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %solution in term of velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(0,yy,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = reshape(InOut,radial,axial);
    %InOut = (InOut>4)/visc + (InOut<4);
    InOut = (InOut>pi)/visc + (InOut<pi);% + (InOut>=3.14 | InOut<=3.15)*2/(1+visc);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    [PX,PY,PI1,PI2] = computeP_AX_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    
    %compute underlying extensional flow
    [ux_ext,uy_ext] = extens_flow(X0',Y0',capillary);
    
    %compute the integral to obtain the velocities
    U = ux_ext - 1/8/pi*GXX*df_x' - 1/8/pi*GXY*df_y' + (1-visc)/8/pi*A11*ux_sol + (1-visc)/8/pi*A12*uy_sol;
    V = uy_ext - 1/8/pi*GYX*df_x' - 1/8/pi*GYY*df_y' + (1-visc)/8/pi*A21*ux_sol + (1-visc)/8/pi*A22*uy_sol;
    
    %compute pressure
    if visc==1
        p = -(PX*df_x'+PY*df_y')/8/pi;
    else
        p = -(PX*df_x'+PY*df_y')/8/pi + (1-visc)./InOut.*(PI1*ux_sol + PI2*uy_sol)/8/pi;
    end
    
    %reshape in the grid shape
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    p = reshape(p,radial,axial);
    
    %moltiply time 1/lambda the points inside the droplet
    U = InOut.*U;
    V = InOut.*V;
    
%     if PlotHere==1
%     
%     figure
%     quiver(X,Y,U,V)
%     axis equal
%     hold on
%     quiver(X,-Y,U,-V,'b')
%     plot(a,b,'k-',a,-b,'k-','LineWidth',2)
%     hold off
%     xlabel('axial direction')
%     ylabel('radial direction')
%     
% %     figure
% %     quiver(Y,-X,V,-U)
% %     axis equal
% %     hold on
% %     quiver(-Y,-X,-V,-U,'b')
% %     plot(b,-a,'r',-b,-a,'r')
% %     hold off
% %     ylabel('axial direction')
% %     xlabel('radial direction')
%     
%     figure
%     contourf(X,Y,abs(U),20)
%     hold on
%     contour(X,-Y,abs(U),20)
%     plot(a,b,'k-',a,-b,'k-','LineWidth',2)
%     axis equal
%     hold off
%     xlabel('axial direction')
%     ylabel('radial direction')
%     title('axial velocity')
%     
%     figure
%     contourf(X,Y,abs(V),20)
%     hold on
%     plot(a,b,'k-',a,-b,'k-','LineWidth',2)
%     contour(X,-Y,abs(V),20)
%     axis equal
%     hold off
%     xlabel('axial direction')
%     ylabel('radial direction')
%     title('radial velocity')
%     
%     end
    
%     figure
%     contour(X,Y,sqrt(U.*U+V.*V),20)
%     hold on
%     plot(a,b,'k-',a,-b,'k-','LineWidth',2)
%     contour(X,-Y,sqrt(U.*U+V.*V),20)
%     axis equal
%     hold off
%     xlabel('axial direction')
%     ylabel('radial direction')
%     title('velocity magnitude')
    
end