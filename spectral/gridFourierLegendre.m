%create grid theta phi for spherical harmonics

function [dirs,w,thetaLine,phiLine] = gridFourierLegendre(nLegendre,nFourier)

%derivation matrices and integration weight
%[theta,D1t,D2t,wt,~,PARAM.PPP] = LegendreIntDiff(0,pi,nLegendre);
%wLegendre = wt;
%PARAM.wT = wLegendre;
%[phi,D1s] = fourdif(nFourier,1);
%[~,D2s] = fourdif(nS,2);
%wS = ws;
%PARAM.wS = ws;

[theta,~,~,wTheta] = LegendreIntDiff(0,pi,nLegendre);
phi = fourdif(nFourier,1);
wPhi = fourierIntWeight(nFourier)';

thetaLine = theta;
phiLine = phi;
[theta,phi] = meshgrid(theta,phi);
theta = theta(:);
phi = phi(:);

dirs = [phi theta];

[wTheta,wPhi] =  meshgrid(wTheta,wPhi);
wTheta = wTheta(:);     wPhi = wPhi(:);

w = wTheta.*wPhi;