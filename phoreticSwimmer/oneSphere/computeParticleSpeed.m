%handle function for one particle

function U = computeParticleSpeed(fMode,PARAM)

%physical grid
theta = linspace(0,pi,PARAM.n+1)';

%detrmine first mode such that the volume is conserved
fVolume = @(unk) ModifyVolumeModes2(theta,fMode,unk,PARAM.V0);
options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
f0 = fsolve(fVolume,1,options);

%compute shape
fMode = [f0; fMode];
r = LegendreBuildShape(theta,fMode,0);
xHere = r'.*cos(theta');
yHere = r'.*sin(theta');
yHere([1 end]) = [0 0];
x{1} = xHere;
y{1} = yHere;

%solve laplace equation
conc = BEM_Laplace(x,y,PARAM);

%compute slip
vSlip = computeSlipVelPhoretic(x{1},y{1},conc(1:PARAM.n),PARAM,1,PARAM.D1{1});

%function profile for BCs
PARAM.velBC{1} = vSlip;

%solve Stokes equation
yStokes = BEM_Stokes(x,y,PARAM);

U = yStokes(end);