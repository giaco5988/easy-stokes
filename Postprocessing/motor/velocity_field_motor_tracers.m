%droplet VELOCITY FIELD visualization of one droplet in a channel (Etienne Lac like)

function [U,V,U0] = velocity_field_motor_tracers(a,b,solution,X,Y,ite,PARAM,LAB,PlotNow,MASK)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
  
  %option
  vort = 0; %plot vorticity
  %LAB = 1;  %stay in motor frame or lab frame or drop frame
  
  m = PARAM.m;    PARAM.q = numel(a)-2-m;   q = PARAM.q;
  fixed_elem = m;

  %compute the spline coeff
  [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(a(m+2:end),b(m+2:end));

  visc = PARAM.visc;
  visc2 = 1;
  visc1 = PARAM.visc*visc2;
  %Ca = PARAM.Ca;
  dt = PARAM.deltaT;
  checkpoint = PARAM.checkpoint;
  
  %compute singularities coordinates
  Xsing = [(a(1:fixed_elem)+a(2:fixed_elem+1))/2 a(fixed_elem+2:end)];
  Ysing = [(b(1:fixed_elem)+b(2:fixed_elem+1))/2 b(fixed_elem+2:end)];
  
  %normla to wall
  den = sqrt(diff(a(1:m+1)).^2+diff(b(1:m+1)).^2);
  r = [diff(b(1:m+1))./den; -diff(a(1:m+1))./den];
    
  %compute gamma with the CURRENT DEFINITION OF Ca
  radius = PARAM.R; %motor average radius
  U = PARAM.mass/4/pi/radius^2;
  
  %surface tension based on mass injection
  gamma = abs(U)*PARAM.visc2/PARAM.Ca;
      
  if PARAM.continuity==0
      PARAM.Ca = 1;
      gamma = PARAM.visc2*PARAM.R/PARAM.Ca;
  end
  
  %BC stresses drop
  [df_x,df_y] = df_surf_tens_spline_symm(a(m+2:end),b(m+2:end),gamma);
  
  %compute normal for nodes
  N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
      -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];
  
  %add special forces
  if PARAM.SpecialForces==1
      
      %compute distance between droplet and wall FOR DROP
      distDWdrop = distDropWall(a(PARAM.m+2:PARAM.q+PARAM.m+2),b(PARAM.m+2:PARAM.q+PARAM.m+2),a(1:PARAM.m+1),b(1:PARAM.m+1));
      
      %add Van der Walls with Hamacker model
      ForceVanDerWalls = -PARAM.strength*distDWdrop.^(-3);
      df_x = df_x + ForceVanDerWalls.*N(1,:);
      df_y = df_y + ForceVanDerWalls.*N(2,:);
      
      %compute integral of VanDerWalls in order to balance zero force
      %condition
      temp = ForceVanDerWalls.*N(1,:);   dl = sqrt(diff(a(fixed_elem+2:end)).^2+diff(b(fixed_elem+2:end)).^2);
      FX = pi*sum((temp(1:end-1).*b(fixed_elem+2:end-1)+temp(2:end).*b(fixed_elem+3:end)).*dl);
      
  elseif PARAM.SpecialForces==2
      
      %compute distance between droplet and wall FOR DROP
      distDWdrop = distDropWall(a(PARAM.m+2:PARAM.q+PARAM.m+2),b(PARAM.m+2:PARAM.q+PARAM.m+2),a(1:PARAM.m+1),b(1:PARAM.m+1));
      
      %compute disjoining pressure
      PI = disj_pressure(distDWdrop,PARAM.cut,PARAM.strength);
      df_x = df_x + PI.*N(1,:);
      df_y = df_y + PI.*N(2,:);
      
      %compute integral of VanDerWalls in order to balance zero force
      %condition
      temp = PI.*N(1,:);   dl = sqrt(diff(a(fixed_elem+2:end)).^2+diff(b(fixed_elem+2:end)).^2);
      FX = pi*sum((temp(1:end-1).*b(fixed_elem+2:end-1)+temp(2:end).*b(fixed_elem+3:end)).*dl);
      
      
  end
  
  %modify BC for BEM
  uxPotential = 0;
  uyPotential = 0;
  fxPotential = 0;
  fyPotential = 0;
  NXX = 0;
  NXY = 0;
  NYX = 0;
  NYY = 0;
  if PARAM.continuity==1
      
      PARAM.x_inj(1) = center_mass(a(m+2:end),b(m+2:end));
      
      for l = 1:numel(PARAM.mass)
  
        %add velocity field and stresses with potential flow
        [ux_temp,uy_temp,nxx_temp,nxy_temp,nyx_temp,nyy_temp] = gf2D_ax_laplace_fs(Xsing,Ysing,PARAM.x_inj(l),PARAM.mass(l));
          
              uxPotential = uxPotential + ux_temp;
              uyPotential = uyPotential + uy_temp;
%               NXX = NXX + PARAM.visc2*[nxx_temp(1:fixed_elem) zeros(1,q+1)];
%               NXY = NXY + PARAM.visc2*[nxy_temp(1:fixed_elem) zeros(1,q+1)];
%               NYX = NYX + PARAM.visc2*[nyx_temp(1:fixed_elem) zeros(1,q+1)];
%               NYY = NYY + PARAM.visc2*[nyy_temp(1:fixed_elem) zeros(1,q+1)];
              NXX = NXX + zeros(1,m+q+1);
              NXY = NXY + zeros(1,m+q+1);
              NYX = NYX + zeros(1,m+q+1);
              NYY = NYY + zeros(1,m+q+1);
              fxPotential = NXX.*[r(1,:) N(1,:)]+NXY.*[r(2,:) N(2,:)];
              fyPotential = NYX.*[r(1,:) N(1,:)]+NYY.*[r(2,:) N(2,:)];
    
      end
      
  end

  %get solution
  solution = solution(1:2*(m+q+1)+1,ite);
  
  %compute splines coordinates
  t = 0:0.1:0.9;
    
  %INTERFACE
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

  %mesh finess per unit length
  %radial = round(resy*yy);
  %axial = round(resx*abs(xx1-xx2));
  [radial,axial] = size(X);
      
  %velocities
  ux = [zeros(m,1)*solution(end); solution(2*m+1:2:2*(m+q)+1)] - uxPotential';
  uy = [zeros(m,1); solution(2*m+2:2:2*(m+q)+2)] - uyPotential';
  
  %solution
  fx = [solution(1:2:2*m-1); df_x'] - fxPotential';
  fy = [solution(2:2:2*m); df_y'] - fyPotential';
    
  X0 = X(:);
  Y0 = Y(:);
  
  numGRID = numel(X0);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Green's functions
  [GXXwall,GXYwall,GYXwall,GYYwall,TXXXwall,TXXYwall,TXYYwall,TYXXwall,TYXYwall,TYYYwall] =...
            Stokes2DAxisIntConst(a(1:m+1),b(1:m+1),X0',Y0');
        
  R1 = repmat(r(1,:),numGRID,1);
  R2 = repmat(r(2,:),numGRID,1);
        
  GXXwall = reshape(GXXwall,numGRID,fixed_elem);
  GXYwall = reshape(GXYwall,numGRID,fixed_elem);
  GYXwall = reshape(GYXwall,numGRID,fixed_elem);
  GYYwall = reshape(GYYwall,numGRID,fixed_elem);
  TXXXwall = reshape(TXXXwall,numGRID,fixed_elem);
  TXXYwall = reshape(TXXYwall,numGRID,fixed_elem);
  TXYYwall = reshape(TXYYwall,numGRID,fixed_elem);
  TYXXwall = reshape(TYXXwall,numGRID,fixed_elem);
  TYXYwall = reshape(TYXYwall,numGRID,fixed_elem);
  TYYYwall = reshape(TYYYwall,numGRID,fixed_elem);
       
  A11wall = TXXXwall.*R1 + TXXYwall.*R2;
  A12wall = TXXYwall.*R1 + TXYYwall.*R2;
  A21wall = TYXXwall.*R1 + TYXYwall.*R2;
  A22wall = TYXYwall.*R1 + TYYYwall.*R2;
        
  [GXXint,GXYint,GYXint,GYYint,TXXXint,TXXYint,TXYXint,TXYYint,TYXXint,TYXYint,TYYXint,TYYYint] =...
       Stokes2DAxisSPlinesLinear(ax,bx,cx,dx,ay,by,cy,dy,X0',Y0');
        
  GXXint = reshape(GXXint,numGRID,q+1);
  GXYint = reshape(GXYint,numGRID,q+1);
  GYXint = reshape(GYXint,numGRID,q+1);
  GYYint = reshape(GYYint,numGRID,q+1);
  TXXXint = reshape(TXXXint,numGRID,q+1);
  TXXYint = reshape(TXXYint,numGRID,q+1);
  TXYXint = reshape(TXYXint,numGRID,q+1);
  TXYYint = reshape(TXYYint,numGRID,q+1);
  TYXXint = reshape(TYXXint,numGRID,q+1);
  TYXYint = reshape(TYXYint,numGRID,q+1);
  TYYXint = reshape(TYYXint,numGRID,q+1);
  TYYYint = reshape(TYYYint,numGRID,q+1);
        
  A11int = TXXXint + TXXYint;
  A12int = TXYXint + TXYYint;
  A21int = TYXXint + TYXYint;
  A22int = TYYXint + TYYYint;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %BCwall
  uxWall = ux(1:fixed_elem);    uyWall = uy(1:fixed_elem);
  fxWall = fx(1:fixed_elem);    fyWall = fy(1:fixed_elem);

  %BC INTERFACE
  uxInter = ux(fixed_elem+1:end);    uyInter = uy(fixed_elem+1:end);
  fxInter = fx(fixed_elem+1:end);    fyInter = fy(fixed_elem+1:end);

  %compute the integral to obtain center of mass velocity
  %compute normal velocity
  %N = [by./sqrt(bx.*bx+by.*by) by(end)+2*cy(end)+3*dy(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      %-bx./sqrt(bx.*bx+by.*by) -bx(end)-2*cx(end)-3*dx(end)/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];
  uNormal = N(1,:)'.*uxInter + N(2,:)'.*uyInter;
  vCentermass = DropVelocityAxis(a(m+2:end),b(m+2:end),uNormal);

  %compute the integral ON WALLS to obtain the velocities
  U = -GXXwall*fxWall-GXYwall*fyWall + visc2*(A11wall*uxWall+A12wall*uyWall);
  V = -GYXwall*fxWall-GYYwall*fyWall + visc2*(A21wall*uxWall+A22wall*uyWall);
  %U = -GXXwall*fxWall-GXYwall*fyWall;
  %V = -GYXwall*fxWall-GYYwall*fyWall;
    
  %compute the integral ON INTERFACE to obtain the velocities
  U = (U-GXXint*fxInter-GXYint*fyInter + (visc2-visc1)*(A11int*uxInter+A12int*uyInter))/8/pi;
  V = (V-GYXint*fxInter-GYYint*fyInter + (visc2-visc1)*(A21int*uxInter+A22int*uyInter))/8/pi;

  %figure out if I'm inside or outside the droplet
  InOut = FigInOut(a(fixed_elem+2:end),b(fixed_elem+2:end),X0,Y0);
  if MASK==1
    InOut = (InOut>pi)*0 + (InOut<pi);
  elseif MASK==0
    InOut = (InOut>pi)/PARAM.visc + (InOut<pi);
  end

  %check if evaluation point is very close to the interface and change it if
  %it is so
  p = 1;
  [X0,Y0,U,V,~,InOut] = close_interface_visu(a(fixed_elem+2:end),b(fixed_elem+2:end),X0,Y0,U,V,p,InOut,solution(2*fixed_elem+1:end),Inf);

  %add flow due to mass injection
  uxField = 0;
  uyField = 0;
  if PARAM.continuity==1
      
      for l = 1:numel(PARAM.mass)
  
        %add velocity field and stresses with potential flow
        [ux_temp,uy_temp] = gf2D_ax_laplace_fs(X0,Y0,PARAM.x_inj(l),PARAM.mass(l));
          
         uxField = uxField + ux_temp;
         uyField = uyField + uy_temp;
    
      end
    
  end
  
  %reshape in the grid shape
  X = reshape(X0,radial,axial);
  Y = reshape(Y0,radial,axial);
  U = reshape(U,radial,axial);
  V = reshape(V,radial,axial);
  if PARAM.continuity==1
    uxField = reshape(uxField,radial,axial);
    uyField = reshape(uyField,radial,axial);
  end
  %p = reshape(p,radial,axial);

  %moltiply time 1/lambda the points inside the droplet
  InOut = reshape(InOut,radial,axial);
  U = InOut.*U + uxField;
  V = InOut.*V + uyField;

  if vort==1
      %compute vorticity with finite differences
      omega = vorticity(axial,radial,U,V,X,Y);
  end

  %axial velocity of the interface (GEOMETRICAL CENTER OF MASS)
  if LAB==0
      if PARAM.lab==1
        U0 = solution(end);
      elseif PARAM.lab==0
        U0 = 0;
      end
  elseif LAB==1
      if PARAM.lab==1
        U0 = 0;
      elseif PARAM.lab==0
        U0 = solution(end);
      end
  elseif LAB==2
      %if PARAM.lab==1
        U0 = vCentermass;
      %elseif PARAM.lab==0
       % U0 = solution(end);
      %end
  end
  %U0 = v_centermass;
  
  if MASK==1
      
      U = InOut.* (U + uxField);
      V = InOut.* (V + uyField);
      U0 = InOut*U0;
  
  end

  if PlotNow==1
      
  figure
  quiver(X,Y,U-U0,V,'b')
  axis equal
  hold on
  plot(xxx,yyy,'k')
  plot(xxx,-yyy,'k')
  plot(a(1:m+1),b(1:m+1),'k')
  plot(a(1:m+1),-b(1:m+1),'k')
  %quiver(Xsing,Ysing,ux',uy')
  %quiver((a(1:end-1)+a(2:end))/2,-(b(1:end-1)+b(2:end))/2,ux',-uy','g')
  quiver(X,-Y,U-U0,-V,'b')
  %plot(PARAM.x_inj,-PARAM.y_inj,'or')
  hold off
  xlabel('axial direction')
  ylabel('radial direction')
  title(['flow field, \lambda=' num2str(visc)])
  
  figure
  streamslice(X,Y,U-U0,V,'b')
  axis equal
  hold on
  plot(xxx,yyy,'k')
  plot(xxx,-yyy,'k')
  plot(a(1:m+1),b(1:m+1),'k')
  plot(a(1:m+1),-b(1:m+1),'k')
  %quiver(Xsing,Ysing,ux',uy')
  %quiver((a(1:end-1)+a(2:end))/2,-(b(1:end-1)+b(2:end))/2,ux',-uy','g')
  quiver(X,-Y,U-U0,-V,'b')
  %normVel = sqrt((U-U0).^2+V.^2);
  %quiver(X,-Y,(U-U0)./normVel,-V./normVel,'b')
  %plot(PARAM.x_inj,-PARAM.y_inj,'or')
  hold off
  xlabel('axial direction')
  ylabel('radial direction')
  title(['flow field, \lambda=' num2str(visc)])

  if vort==1

    figure
    contourf(X,Y,omega,200,'LineStyle','none')
    axis equal
    %axis([-yy yy xx1 xx2])
    hold on
    %plot(Y0,X0,'or')
    contourf(X,-Y,omega,200,'LineStyle','none')
    plot(xxx,yyy,'k',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
    plot(a(1:m+1),b(1:m+1),'k')
    plot(a(1:m+1),-b(1:m+1),'k')
    hold off
    colorbar
    ylabel('axial direction')
    xlabel('radial direction')
    title(['Vorticity at t=' num2str(ite*dt*checkpoint) ', \lambda=' num2str(visc) ', Ca=' num2str(Ca)])

  end
  
  end
    
end