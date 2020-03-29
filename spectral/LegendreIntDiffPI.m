%compute Legendre for integration and derivation for n point grid

function [t,D1,D2,w,manyWG,PPP,PPPbc] = LegendreIntDiffPI(x0,x1,n)

%differentiation matrix and points for Legendre
L = x1-x0;
scale = 2/L + x0;
[x,w] = legs(n);
w = w';
t=(x+1)/scale;
D1=legsdiff(n,x);
D1 = D1*scale;
D2 = D1*D1;

%integrationn weights and points for Legendre
w = w'/scale;
WG = w;
prepareWG = repmat(WG',2,1);
manyWG = repmat(prepareWG(:)',2*n,1);

%all legendre polynomia that I need
PPP = zeros(n,n);
PPPbc = zeros(n,2);
for i = 1:n
    P = legendre(i-1,2*t/L-1);
    P = P(1,:)';   %legendre polynomia
    PPP(i,:) = P;
    
    %for boundary conditions
    Pbc = legendre(i-1,2*[0; 1]-1);
    PPPbc(i,:) = Pbc(1,:)';
end