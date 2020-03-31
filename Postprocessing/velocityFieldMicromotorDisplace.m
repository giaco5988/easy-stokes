%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,U0,p] = velocityFieldMicromotorDisplace(aDrop,bDrop,aWall,bWall,solution,visc,X,Y,PARAM,LAB,MASK,decompose)
  
  if PARAM.resizeDFXimplicit==1
      
     dfStar = solution(end);
     solution = solution(1:end-1);
      
  end
  
  Xsing = [(aWall(1:end-1)+aWall(2:end))/2 aDrop];
  Ysing = [(bWall(1:end-1)+bWall(2:end))/2 bDrop];

  %often used variables
  m = numel(aWall)-1;  q = numel(aDrop)-1;
  fixed_elem = m;
  
  %mass injection
  radius = PARAM.R;
  PARAM.Qsource = PARAM.Ca*radius^2*1/PARAM.visc2;
  
  %viscosities
  visc2 = 1;
  visc1 = visc;
  
  %normal vector
  noWall = sqrt(diff(bWall).^2+diff(aWall).^2); noDrop = sqrt(diff(bDrop).^2+diff(aDrop).^2);
  r = [-diff(bWall)./noWall -diff(bDrop)./noDrop; diff(aWall)./noWall diff(aDrop)./noDrop];

  %mesh finess per unit length
  [radial,axial] = size(X);
  
  %compute stresses on droplet form geometry
  %compute spline coefficient and curvature
  [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric(aDrop, bDrop);
  
  %compute normal for nodes
  N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
        -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];
    
  %compute curveture
  K1 = curv_spline2(bx,by,cx,cy,dx,dy);
  K2 = N(2,:)./bDrop;
  K2([1 end]) = K1([1 end]);
  K = K1+K2;
  df = K;
  
  % add disjoinig pressure
  if PARAM.repulsive==2;
      %compute distance wall-drop
         distDrop = distWallDrop(aWall,bWall,aDrop,bDrop);
              
         if min(distDrop)<PARAM.repulsiveOn
            dfRepulsive = repulsiveStress(PARAM,distDrop);
            
            df = df + dfRepulsive;
         end
         
         %stresses from the BCs
          dfXbc = df.*N(1,:);
          dfYbc = df.*N(2,:);
     
  elseif PARAM.repulsive==3
          
        %compute distance wall-drop
        distDrop = distWallDrop(aWall,bWall,aDrop,bDrop);
          
        %display('Hamacker')
          
        dfRepulsive = disjoiningPressure(PARAM,distDrop);
        df = df + dfRepulsive;
        
        %stresses from the BCs
        dfXbc = df.*N(1,:);
        dfYbc = df.*N(2,:);
        
  elseif PARAM.repulsive==4
      
        %compute distance wall-drop
        distDrop = distWallDrop(aWall,bWall,aDrop,bDrop);
          
        %display('Hamacker')
          
        dfRepulsive = disjoiningPressure(PARAM,distDrop);
        
        %stresses from the BCs
        dfXbc = df.*N(1,:);
        dfYbc = (df+dfRepulsive).*N(2,:);
        
  elseif PARAM.repulsive==5
          
          %compute distance wall-drop
          distDrop = distWallDrop(aWall,bWall,aDrop,bDrop);
          
          display('Double layer')
          
          dfRepulsive = disjoiningPressureDoubleLayer(PARAM,distDrop);
          
          df = df + dfRepulsive;
          
          dfXbc = df.*N(1,:);
          dfYbc = df.*N(2,:);
          
  elseif PARAM.repulsive==6
          
          display('Double layer direction')
          
          [dfXtoDrop,dfYtoDrop,dfXtoWall,dfYtoWall] = disjoiningPressureDoubleLayerDirection(Xsing(1:m),Ysing(1:m),aDrop,bDrop,PARAM);
          
          dfXbc = df.*N(1,:) + dfXtoDrop';
          dfYbc = df.*N(2,:) + dfYtoDrop';
        
  elseif PARAM.repulsive==0
      
          dfXbc = df.*N(1,:);
          dfYbc = df.*N(2,:);
         
  end
  
  if PARAM.resizeDFX==1
          
          fStressX = @(unk) rescaleFunction(aDrop,bDrop,dfXbc',unk,0);
          options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
          deltaFX = fsolve(fStressX,1,options);
          
          dfXbc = dfXbc + deltaFX;
          
          %DropForce = int_axis_spline_symmetric(aDrop,bDrop,dfXbc');
          %DropForce = forceOnDrop(aDrop',bDrop',K',1);
          %display(DropForce)
          
  elseif PARAM.resizeDFX==2
          
          %compute compensation normal
          fStressX = @(unk) rescaleFunction2(aDrop,bDrop,dfXbc,unk,0,N(1,:));
          options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
          deltaFX = fsolve(fStressX,1,options);
          
          %add stress components
          dfXbc = dfXbc + deltaFX.*N(1,:).^2;
          
  end
  
  if PARAM.resizeDFXimplicit==1
      
      dfXbc = dfXbc + dfStar;
      
  end
  
  %velocities from the BCs
  uXbc = zeros(m,1);
  uYbc = zeros(m,1);
  
  %stresses from the solution
  dfX = solution(1:2:2*m-1);
  dfY = solution(2:2:2*m);
  
  if PARAM.MotorFree==3
      
      dfX = dfX + dfXtoWall;
      dfY = dfY + dfYtoWall;
      
  end
  
  %velocities from the solution
  uX = solution(2*m+1:2:end-1);
  uY = solution(2*m+2:2:end);
  
  if decompose==1
      %only drop
      dfX = zeros(m,1);
      dfY = zeros(m,1);
      PARAM.Qsource = 0;
  elseif decompose==2
      %only motor
      dfXbc = zeros(1+q,1)';
      dfYbc = zeros(1+q,1)';
      uX = zeros(1+q,1);
      uY = zeros(1+q,1);
      PARAM.Qsource = 0;
  elseif decompose==3
      %only mass
      dfX = zeros(m,1);
      dfY = zeros(m,1);
      dfXbc = zeros(1+q,1)';
      dfYbc = zeros(1+q,1)';
      uX = zeros(1+q,1);
      uY = zeros(1+q,1);
  elseif decompose==4
      %bubble and motor
      PARAM.Qsource = 0;
  elseif decompose==5
      %bubble and mass
      dfX = zeros(m,1);
      dfY = zeros(m,1);
  elseif decompose==6
      %motor and mass
      dfXbc = zeros(1+q,1)';
      dfYbc = zeros(1+q,1)';
      uX = zeros(1+q,1);
      uY = zeros(1+q,1);
  end
  
  %create all values onthe bounadries
  dfx = [dfX; dfXbc'];
  dfy = [dfY; dfYbc'];
  ux = [uXbc; uX];
  uy = [uYbc; uY];
    
  X0 = X(:);
  Y0 = Y(:);
  
  numField = numel(X0);

  %compute the necessary Green's function
  [GXXwall,GXYwall,GYXwall,GYYwall,TXXXwall,TXXYwall,TXYYwall,TYXXwall,TYXYwall,TYYYwall] =...
            Stokes2DAxisIntConst(aWall,bWall,X0',Y0');
        
  R1 = repmat(r(1,1:fixed_elem),numField,1);
  R2 = repmat(r(2,1:fixed_elem),numField,1);
        
  GXXwall = reshape(GXXwall,numField,fixed_elem);
  GXYwall = reshape(GXYwall,numField,fixed_elem);
  GYXwall = reshape(GYXwall,numField,fixed_elem);
  GYYwall = reshape(GYYwall,numField,fixed_elem);
  TXXXwall = reshape(TXXXwall,numField,fixed_elem);
  TXXYwall = reshape(TXXYwall,numField,fixed_elem);
  TXYYwall = reshape(TXYYwall,numField,fixed_elem);
  TYXXwall = reshape(TYXXwall,numField,fixed_elem);
  TYXYwall = reshape(TYXYwall,numField,fixed_elem);
  TYYYwall = reshape(TYYYwall,numField,fixed_elem);
        
  A11wall = TXXXwall.*R1 + TXXYwall.*R2;
  A12wall = TXXYwall.*R1 + TXYYwall.*R2;
  A21wall = TYXXwall.*R1 + TYXYwall.*R2;
  A22wall = TYXYwall.*R1 + TYYYwall.*R2;
        
  [GXXint,GXYint,GYXint,GYYint,TXXXint,TXXYint,TXYXint,TXYYint,TYXXint,TYXYint,TYYXint,TYYYint] =...
            Stokes2DAxisIntLinearProva(aDrop,bDrop,X0',Y0');
        
  GXXint = reshape(GXXint,numField,q+1);
  GXYint = reshape(GXYint,numField,q+1);
  GYXint = reshape(GYXint,numField,q+1);
  GYYint = reshape(GYYint,numField,q+1);
  TXXXint = reshape(TXXXint,numField,q+1);
  TXXYint = reshape(TXXYint,numField,q+1);
  TXYXint = reshape(TXYXint,numField,q+1);
  TXYYint = reshape(TXYYint,numField,q+1);
  TYXXint = reshape(TYXXint,numField,q+1);
  TYXYint = reshape(TYXYint,numField,q+1);
  TYYXint = reshape(TYYXint,numField,q+1);
  TYYYint = reshape(TYYYint,numField,q+1);
        
  A11int = TXXXint + TXXYint;
  A12int = TXYXint + TXYYint;
  A21int = TYXXint + TYXYint;
  A22int = TYYXint + TYYYint;
        
  GXX = [GXXwall GXXint];
  GXY = [GXYwall GXYint];
  GYX = [GYXwall GYXint];
  GYY = [GYYwall GYYint];
  A11 = [A11wall A11int];
  A12 = [A12wall A12int];
  A21 = [A21wall A21int];
  A22 = [A22wall A22int];
  %[PX,PY] = computeP_visu_const_StokesAX(a,b,X0,Y0);
  
  %compute underlying mass injection
  [uMass,vMass] = massInjection(X0-PARAM.Xinj,Y0,PARAM.Qsource);
  %uMass = 0;
  %vMass = 0;
  
  %compute the integral to obtain the velocities
  singleX = -(GXX*dfx + GXY*dfy);   singleY = -(GYX*dfx + GYY*dfy);
  doubleX = (visc2-visc1)*(A11(:,m+1:end)*ux(m+1:end) + A12(:,m+1:end)*uy(m+1:end));
  doubleY = (visc2-visc1)*(A21(:,m+1:end)*ux(m+1:end) + A22(:,m+1:end)*uy(m+1:end));
  U = uMass*visc1 + (singleX + doubleX)/8/pi;
  V = vMass*visc1 + (singleY + doubleY)/8/pi;
    
  %compute pressure
%   if visc==1
%       p = -(PX*dfX+PY*dfY)/8/pi;
%   else
%       p = -(PX*dfX+PY*dfY)/8/pi;
%   end
    
  %reshape in the grid shape
  U = reshape(U,radial,axial);
  V = reshape(V,radial,axial);
  %p = reshape(p,radial,axial);
  
  %set to zero velocities in wall and apply viscosity ratio
  InOutDrop = FigInOut(aDrop,bDrop,X0,Y0);
  InOutWall = FigInOut(flip(aWall),flip(bWall),X0,Y0);
    
  InOutDrop = reshape(InOutDrop,radial,axial);
  InOutWall = reshape(InOutWall,radial,axial);
    
  U = (InOutDrop<pi).*U + (InOutDrop>pi).*U/PARAM.visc;
  V = (InOutDrop<pi).*V + (InOutDrop>pi).*V/PARAM.visc;
  U = (InOutWall<pi).*U;    V = (InOutWall<pi).*V;
  
  %set to zero velocities in wall and in drop
  if MASK==1
 
    U = (InOutDrop<pi&InOutWall<pi).*U;
    V = (InOutDrop<pi&InOutWall<pi).*V;
    
  end
    
  %moltiply time 1/lambda the points inside the droplet
%   U = InOut.*U;
%   V = InOut.*V;
  
  %choose frame of reference
  if LAB==0
      
       U0 = solution(end);
      
  elseif LAB==1
      
       U0 = 0;
 
  elseif LAB==2
      
       uNormal = N(1,:)'.*uX + N(2,:)'.*uY;
       vCentermass = DropVelocityAxis(aDrop,bDrop,uNormal);
  
       U0 = vCentermass;
  end
  
end