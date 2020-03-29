%compute curvature in axisymmetric

function [K1,K2,K] = curvatureAxisSpectralXY(x,y,PARAM)

%compute rho in symmetry axis
D1 = PARAM.D1;
D2 = PARAM.D2;

%compute geomtrical derivaties
xp = D1*x;    yp = D1*y;
xpp = D2*x;    ypp = D2*y;

%compute normal vector
h = (xp.^2+yp.^2).^(0.5);
ny = -xp./h;

%compute curvature in meridional plane
K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);

%compute curvature in aximuthal direction
K2 = ny./y;

%sum of the two curavtures
K = K1+K2;