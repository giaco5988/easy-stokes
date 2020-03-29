%derivatives on a 2D grid

function [xt,xs,yt,ys,zt,zs] = firstDerivative2Dgrid(x,y,z,PARAM)

%differentiation matrices
nT = PARAM.nT;
nS = PARAM.nS;
D1t = PARAM.D1t;
D1s = PARAM.D1s;

%compute derivative x
x = reshape(x,nS,nT);
xt = (D1t*x')';
xt = xt(:);
xs = (D1s*x);
xs = xs(:);
xt = reshape(xt,nS,nT);
xt = xt(:);
    
%compute derivative y
y = reshape(y,nS,nT);
yt = (D1t*y')';
yt = yt(:);
ys = (D1s*y);
ys = ys(:);
yt = reshape(yt,nS,nT);
yt = yt(:);
    
%compute derivative z
z = reshape(z,nS,nT);
zt = (D1t*z')';
zt = zt(:);
zs = (D1s*z);
zs = zs(:);
zt = reshape(zt,nS,nT);
zt = zt(:);
    
    
    
    
    
    
    
    
    
    
    
    
    
    