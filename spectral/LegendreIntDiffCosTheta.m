
%compute Legendre for integration and derivation for n point grid

function [cosTheta,D1,D2,w,manyWG,PPP,PPPbc] = LegendreIntDiffCosTheta(n)

%differentiation matrix and points for Legendre
L = 2;
scale = -2/L;
[cosTheta,w] = legs(n);
cosTheta = cosTheta/scale;
w = w'/scale;
D1=legsdiff(n,cosTheta);
D1 = D1*scale;
D2 = D1*D1;
theta = acos(cosTheta);
w = -w./sin(theta');
D2 = (D1.*repmat(-cos(theta),1,n)+D2.*repmat(sin(theta).^2,1,n));
D1 = D1.*repmat(-sin(theta),1,n);

%integrationn weights and points for Legendre
WG = w;
prepareWG = repmat(WG',2,1);
manyWG = repmat(prepareWG(:)',2*n,1);

%all legendre polynomia that I need
PPP = zeros(n,n);
PPPbc = zeros(n,2);
for i = 1:n
    P = legendre(i-1,cosTheta);
    P = P(1,:)';   %legendre polynomia
    PPP(i,:) = P;
end