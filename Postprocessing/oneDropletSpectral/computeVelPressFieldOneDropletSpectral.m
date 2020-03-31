%compute velocity field for one droplet

function [u,v,pressure,K] = computeVelPressFieldOneDropletSpectral(x,y,x0,y0,PARAM)

%derivatives
D1 = PARAM.D1;

[row,col] = size(x0);
x0 = x0(:);
y0 = y0(:);
lambda = PARAM.visc;

%compute green functions
X0 = repmat(x0,1,numel(x));   Y0 = repmat(y0,1,numel(y));
X = repmat(x',numel(x0),1); Y = repmat(y',numel(y0),1);
[SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] = sgf_ax_fs_vect3 (X,Y,X0,Y0);
[LX,LY,NXX,NXY,NYX,NYY] = kernelPressureAxis(X,Y,X0,Y0);

%prepare for integration
xp = D1*x;  yp = D1*y;
metrics = sqrt(xp.^2+yp.^2);
metrics = repmat(metrics',numel(x0),1);
manyW = repmat(PARAM.WG',numel(x0),1);
manyW = manyW.*metrics;
  
%compute solution and stresses
here = pwd;
cd('/Users/Giacomo/Documents/MATLAB/droplet_simulations/dropSpectral')
[sol,nx,ny,~,~,~,~,K] = bemSpectralXY(x,y,PARAM);
cd(here)
ux = sol(1:2:end-1);    uy = sol(2:2:end);
fx = K.*nx;     fy = K.*ny;

%compute velocity field
u = (-manyW.*SXX*fx - manyW.*SXY*fy + (1-lambda)*(manyW.*QXXX*(ux.*nx) + manyW.*QXXY*(ux.*ny) + manyW.*QXYX*(uy.*nx)  + manyW.*QXYY*(uy.*ny)))/8/pi;
v = (-manyW.*SYX*fx - manyW.*SYY*fy + (1-lambda)*(manyW.*QYXX*(ux.*nx) + manyW.*QYXY*(ux.*ny) + manyW.*QYYX*(uy.*nx)  + manyW.*QYYY*(uy.*ny)))/8/pi;

%add undelying flow
if PARAM.BC==1 || PARAM.BC==3
      %add underlying flow
      if PARAM.BC==1    %linear extensional flow
          [uADD,vADD] = extens_flow(x0,y0,PARAM.Ca);
      elseif PARAM.BC==3    %linear extensional flow
          [uADD,vADD] = nonLinearExtensFlow(x0,y0,PARAM.Ca,PARAM.CaNL);
      end
      u = u + uADD';
      v = v + vADD';
end

%compute pressure on the grid
singleLayer = -(manyW.*LX*fx + manyW.*LY*fy);
doubleLayer = (1-lambda)*(manyW.*NXX*(ux.*nx) + manyW.*NYX*(uy.*nx) + manyW.*NXY*(ux.*ny) + manyW.*NYY*(uy.*ny));
pressure = singleLayer/8/pi + doubleLayer/8/pi;

%return to right format
u = reshape(u,row,col);
v = reshape(v,row,col);
pressure = reshape(pressure,row,col);







