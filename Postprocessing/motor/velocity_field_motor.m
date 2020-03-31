%droplet VELOCITY FIELD visualization of one droplet in a channel (Etienne Lac like)

function [U,V,X,Y] = velocity_field_motor(risa,risb,solution,xx1,xx2,yy,resx,resy,ite,PARAM,LAB)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
  
  %option
  vort = 0; %plot vorticity
  %LAB = 1;  %stay in motor frame or lab frame or drop frame
  
  this_folder = pwd;
  mex_files = '~/Documents/MATLAB/mex_files';
  
  addpath(mex_files);
  cd(mex_files);
  mex Stokes2DAxisSPlinesLinear.c
  mex Stokes2DAxisIntConst.c
  cd(this_folder)

  m = PARAM.m;    q = PARAM.q;
  fixed_elem = m;
  
  a = risa(1:m+q+2,ite)';
  b = risb(1:m+q+2,ite)';
  aAfter = risa(1:m+q+2,ite+1)';
  bAfter = risb(1:m+q+2,ite+1)';
  xcm = center_mass(a(2+m:end),b(2+m:end)); xcmAfter = center_mass(aAfter(2+m:end),bAfter(2+m:end));
  vCentermass = (xcmAfter-xcm)/PARAM.deltaT/PARAM.checkpoint;

  visc = PARAM.visc;
  visc2 = 1;
  visc1 = PARAM.visc*visc2;
  Ca = PARAM.Ca;
  dt = PARAM.deltaT;
  checkpoint = PARAM.checkpoint;

  %a = a(:,ite-1)';
  %b = b(:,ite-1)';
  
  %compute singularities coordinates
  Xsing = [(a(1:fixed_elem)+a(2:fixed_elem+1))/2 a(fixed_elem+2:end)];
  Ysing = [(b(1:fixed_elem)+b(2:fixed_elem+1))/2 b(fixed_elem+2:end)];
  
  %normla to wall
  den = sqrt(diff(a(1:m+1)).^2+diff(b(1:m+1)).^2);
  r = [diff(b(1:m+1))./den; -diff(a(1:m+1))./den];
    
  %compute gamma with the CURRENT DEFINITION OF Ca
  radius = PARAM.alpha*PARAM.R; %droplet radius
  U = PARAM.mass/4/pi/radius^2;
  
  %surface tension based on mass injection
  gamma = abs(U)*PARAM.visc2/PARAM.Ca;
      
  if PARAM.continuity==0
      gamma = PARAM.visc2*PARAM.R/PARAM.Ca;
  end

  solution = solution(1:2*(m+q+1)+1,ite);
  %a_before = a(:,ite-2)';
  %b_before = b(:,ite-2)';
    
  %figure out center of mass velocity
  %xcm_old = center_mass(a_before,b_before);
  %xcm = center_mass(a,b);
  %v_centermass = (xcm-xcm_old)/dt/checkpoint;
  
  %compute the spline coeff
  [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(a(m+2:end), b(m+2:end));
  
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
  radial = round(resy*yy);
  axial = round(resx*abs(xx1-xx2));
  
  %BC stresses drop
  [df_x,df_y] = df_surf_tens_spline_symm(a(m+2:end),b(m+2:end),gamma);

      
  %velocities
  ux = [zeros(m,1); solution(2*m+1:2:2*(m+q)+1)];
  uy = [zeros(m,1); solution(2*m+2:2:2*(m+q)+2)];
    
  %solution
  fx = [solution(1:2:2*m-1); df_x'];
  fy = [solution(2:2:2*m); df_y'];

  %create the grid
  x = linspace(xx1,xx2,axial);
  y = linspace(0,yy,radial);

  [X,Y] = meshgrid(x,y);
    
  X0 = X(:);
  Y0 = Y(:);
  
  numGRID = numel(X0);

  %compute the necessary Green's function
%   [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu(a,b,X0,Y0);
%   %[Px,Py] = computeP(PARAM.x_inj,PARAM.y_inj,X0,Y0);
%   
%   r = [-ones(1,n) zeros(1,m) ones(1,j); zeros(1,n) -ones(1,m) zeros(1,j)];
%     
%   R1 = repmat(r(1,1:end),numel(X0),1);
%   R2 = repmat(r(2,1:end),numel(X0),1);
%     
%   A11 = TXXX.*R1+TXXY.*R2;
%   A12 = TXYX.*R1+TXYY.*R2;
%   A21 = TYXX.*R1+TYXY.*R2;
%   A22 = TYYX.*R1+TYYY.*R2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
%   GXX = [GXXwall GXXint];
%   GXY = [GXYwall GXYint];
%   GYX = [GYXwall GYXint];
%   GYY = [GYYwall GYYint];
%   A11 = [A11wall A11int];
%   A12 = [A12wall A12int];
%   A21 = [A21wall A21int];
%   A22 = [A22wall A22int];
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %BCwall
  uxWall = ux(1:fixed_elem);    uyWall = uy(1:fixed_elem);
  fxWall = fx(1:fixed_elem);    fyWall = fy(1:fixed_elem);

  %BC INTERFACE
  uxInter = ux(fixed_elem+1:end);    uyInter = uy(fixed_elem+1:end);
  fxInter = fx(fixed_elem+1:end);    fyInter = fy(fixed_elem+1:end);

  %compute the integral to obtain the velocities
  %U = 

  %compute the integral ON WALLS to obtain the velocities
  U = -GXXwall*fxWall-GXYwall*fyWall + visc2*(A11wall*uxWall+A12wall*uyWall);
  V = -GYXwall*fxWall-GYYwall*fyWall + visc2*(A21wall*uxWall+A22wall*uyWall);
    
  %compute the integral ON INTERFACE to obtain the velocities
  U = (U-GXXint*fxInter-GXYint*fyInter + (visc2-visc1)*(A11int*uxInter+A12int*uyInter))/8/pi;
  V = (V-GYXint*fxInter-GYYint*fyInter + (visc2-visc1)*(A21int*uxInter+A22int*uyInter))/8/pi;

  %figure out if I'm inside or outside the droplet
  InOut = FigInOut(a(fixed_elem+2:end),b(fixed_elem+2:end),X0,Y0);
  InOut = (InOut>pi)/PARAM.visc + (InOut<pi);

  %check if evaluation point is very close to the interface and change it if
  %it is so
  p = 1;
  [X0,Y0,U,V,~,InOut] = close_interface_visu(a(fixed_elem+2:end),b(fixed_elem+2:end),X0,Y0,U,V,p,InOut,solution(2*fixed_elem+1:end),2);

  %reshape in the grid shape
  X = reshape(X0,radial,axial);
  Y = reshape(Y0,radial,axial);
  U = reshape(U,radial,axial);
  V = reshape(V,radial,axial);
  %p = reshape(p,radial,axial);

  %moltiply time 1/lambda the points inside the droplet
  InOut = reshape(InOut,radial,axial);
  U = InOut.*U;
  V = InOut.*V;

  if vort==1
      %compute vorticity with finite differences
      omega = vorticity(axial,radial,U,V,X,Y);
  end
  
  %velocity of this single point IT DOESN'T GIVE THE VELOCITY OF THE CENTER
  %OF MASS!!!
  %U0 = UV_1point_1drop_channel(a,b,ux,uy,fx,fy,xcm,0,PARAM);

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