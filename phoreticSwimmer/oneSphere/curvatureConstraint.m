%handle function for one particle

function [c,ceq] = curvatureConstraint(fMode,PARAM)

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

%compyte curvature
curv = curv_spline(xHere,yHere);

%disequality constraint
c = max(abs(curv))-PARAM.maxCurv;

%equality constraint
ceq = [];











