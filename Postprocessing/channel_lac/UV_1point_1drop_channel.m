%compute velocity field of few points of droplet in channel

function [U,V] = UV_1point_1drop_channel(a,b,ux,uy,fx,fy,X0,Y0,PARAM)

  %geometry parameters
  n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q= PARAM.q;
  fixed_elem = n+m+j;
  numGRID = numel(X0);
  
  %compute the spline coeff
  [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(a(n+m+j+2:end), b(n+m+j+2:end));
  
  %normla to wall
  r = [-ones(1,n) zeros(1,m) ones(1,j); zeros(1,n) -ones(1,m) zeros(1,j)];

  [GXXwall,GXYwall,GYXwall,GYYwall,TXXXwall,TXXYwall,TXYYwall,TYXXwall,TYXYwall,TYYYwall] =...
            Stokes2DAxisIntConst(a(1:n+m+j+1),b(1:n+m+j+1),X0',Y0');
        
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
        
  GXX = [GXXwall GXXint];
  GXY = [GXYwall GXYint];
  GYX = [GYXwall GYXint];
  GYY = [GYYwall GYYint];
  A11 = [A11wall A11int];
  A12 = [A12wall A12int];
  A21 = [A21wall A21int];
  A22 = [A22wall A22int];
  
  visc1 = PARAM.visc;   visc2 = 1;
  
  %compute the integral to obtain the velocities
  U = -GXX*fx-GXY*fy + (visc2-visc1)*(A11*ux+A12*uy);
  V = -GYX*fx-GYY*fy + (visc2-visc1)*(A21*ux+A22*uy);

end

