%compute Legendre for integration and derivation for n point grid

function [t,D1,D2,w,manyWG,PPP] = LegendreIntDiffTemp(x1,n)

%differentiation matrix and points for Legendre
x0 = -x1;
L = x1-x0;
scale = 2/L;
[x,w] = legs(n);
w = w';
t = x/scale;
D1 = legsdiff(n,x);
D1 = D1*scale;
D2 = D1*D1;

%integrationn weights and points for Legendre
w = w'/scale;
WG = w;
prepareWG = repmat(WG',2,1);
manyWG = repmat(prepareWG(:)',2*n,1);

%all legendre polynomia that I need
PPP = zeros(n,n);
for i = 1:n
    P = legendre(i-1,t*scale);
    P = P(1,:)';   %legendre polynomia
    PPP(i,:) = P;
end