%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p] = velocityFieldNoBubble(a,b,solution,visc,X,Y,PARAM)
  
  %mesh finess per unit length
  [radial,axial] = size(X);
    
  %stresses from the solution
  dfX = solution(1:2:end-1);
  dfY = solution(2:2:end);
    
  X0 = X(:);
  Y0 = Y(:);

  %compute the necessary Green's function
  %[GXX,GXY,GYX,GYY] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
  [GXX,GXY,GYX,GYY] = computeGT_visu(a,b,X0,Y0);
  [PX,PY] = computeP_visu_const_StokesAX(a,b,X0,Y0);
    
  %compute underlying extensional flow
  %PARAM.Qsource = 0;
  [uMass,vMass] = massInjection(X0-PARAM.Xinj,Y0,PARAM.Qsource);
    
  %compute the integral to obtain the velocities
  U = uMass - (GXX*dfX + GXY*dfY)/8/pi;
  V = vMass - (GYX*dfX + GYY*dfY)/8/pi;
  %U = (GXX*dfX + GXY*dfY)/8/pi;
  %V = (GYX*dfX + GYY*dfY)/8/pi;
    
  %compute pressure
  if visc==1
      p = -(PX*dfX+PY*dfY)/8/pi;
  else
      p = -(PX*dfX+PY*dfY)/8/pi;
  end
    
  %reshape in the grid shape
  U = reshape(U,radial,axial);
  V = reshape(V,radial,axial);
  p = reshape(p,radial,axial);
    
  %moltiply time 1/lambda the points inside the droplet
  %U = InOut.*U;
  %V = InOut.*V;
    
end