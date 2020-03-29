%integrate a certain quantity on a line in axysmmetrix

function F = axisIntSpectralCurvilinear(x,y,f,SPECTRAL)

%integration and derivation
w = SPECTRAL.WG';
D1 = SPECTRAL.D1;

%metric term
xp = D1*x;  yp = D1*y;
h = sqrt(xp.^2+yp.^2);

F = 2*pi*w*(y.*h.*f);