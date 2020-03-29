%legendre interpolation

function x = interpLegendreFourierGrid(u,v,fT,an,bn)

n = sqrt(numel(fT));
%m = numel(v);

%build legendre
PPP = lepolym(n-1,cos(u));
xLegendre = LegendreBuild2DgridInterp(fT,PPP);

%build Fourier
xFourier = FourierBuild2DgridInterp(an,bn,v);

x = xFourier.*xLegendre;