%compute Chebyshev for integration and derivation for n point grid

function [t,D1,D2,w,manyWG,TTT] = ChebyshevIntDiff(x0,x1,n)

%differentiation matrix and points Chebyshev-Lobatto
L = x1-x0;
scale = -2/L + x0;
[t,DDD] = chebdif(n,2);
t=(t-1)/scale;
D1 = DDD(:,:,1)*scale;
D2 = DDD(:,:,2)*scale^2;

%integrationn weights and points Clenshaw-Curtis
[~,w] = fclencurt(n,x0,x1);
WG = w;
prepareWG = repmat(WG',2,1);
manyWG = repmat(prepareWG(:)',2*n,1);

%all legendre polynomia that I need
TTT = zeros(n,n);
PPPbc = zeros(n,2);
for i = 1:n
    P = legendre(i-1,2*t/L-1);
    P = P(1,:)';   %legendre polynomia
    TTT(i,:) = P;
    
    %for boundary conditions
    Pbc = legendre(i-1,2*[0; 1]-1);
    PPPbc(i,:) = Pbc(1,:)';
end