%compute force acting on droplet in curvilinear coordinates

function F = ForceOnDropSpectralCurvilinear(x,y,Ca,SPECTRAL)

%integration and derivation
D1 = SPECTRAL.D1;  D2 = SPECTRAL.D2;

%compute geomtrical derivaties
xp = D1*x;    yp = D1*y;
xpp = D2*x;    ypp = D2*y;
  
%compute normal vector
h = (xp.^2+yp.^2).^(0.5);
nx = yp./h;
ny = -xp./h;
  
%compute curvature in meridional plane
K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
  
%compute curvature in aximuthal direction
K2 = ny./y;
K2([1 end]) = K1([1 end]);
  
%total curvature
K = K1+K2;
df = K*Ca;

%compute axial component of stresses
fx = df.*nx;

%compute integral of fx in axisymmetric
F = axisIntSpectralCurvilinear(x,y,fx,SPECTRAL);