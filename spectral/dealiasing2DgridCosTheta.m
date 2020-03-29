%compute modes coefficients from grid points in curvilinear coordinates

function x = dealiasing2DgridCosTheta(theta,x,PARAM)

nS = PARAM.nS;
nT = PARAM.nT;

x = reshape(x,nS,nT);

%compute Legendre Mode
fT = LegendreSerie2DgridCosTheta2(repmat(theta',nT,1),x,PARAM.PPP,PARAM.wT);

%dealiasing Legendre
fT(:,PARAM.dealiasingT+1:end) = 0;
       
%build again
x = LegendreBuild2Dgrid(fT,PARAM.PPP);

%compute Fourier Modes
[an,bn] = FourierSerie2Dgrid(x,PARAM.wS,PARAM.s);

%dealiasing Fourier
an(PARAM.dealiasingS+1:end,:) = 0;
bn(PARAM.dealiasingS+1:end,:) = 0;

%build again Fourier
x = FourierBuild2Dgrid(an,bn,PARAM.s);

x = x(:);





