%derivatives on a 2D grid

function [xt,xs,xtt,xts,xss,yt,ys,ytt,yts,yss,zt,zs,ztt,zts,zss] = derivatives2Dgrid(x,y,z,PARAM)

%differentiation matrices
nT = PARAM.nT;
nS = PARAM.nS;
D1t = PARAM.D1t;
D1s = PARAM.D1s;
D2t = PARAM.D2t;
D2s = PARAM.D2s;

%compute derivative x
x = reshape(x,nS,nT);
xt = (D1t*x')';
xt = xt(:);
xtt = (D2t*x')';
xtt = xtt(:);
xs = (D1s*x);
xs = xs(:);
xss = (D2s*x);
xss = xss(:);
xt = reshape(xt,nS,nT);
xts = (D1s*xt);
xts = xts(:);
xt = xt(:);
    
%compute derivative y
y = reshape(y,nS,nT);
yt = (D1t*y')';
yt = yt(:);
ytt = (D2t*y')';
ytt = ytt(:);
ys = (D1s*y);
ys = ys(:);
yss = (D2s*y);
yss = yss(:);
yt = reshape(yt,nS,nT);
yts = (D1s*yt);
yts = yts(:);
yt = yt(:);
    
%compute derivative z
z = reshape(z,nS,nT);
zt = (D1t*z')';
zt = zt(:);
ztt = (D2t*z')';
ztt = ztt(:);
zs = (D1s*z);
zs = zs(:);
zss = (D2s*z);
zss = zss(:);
zt = reshape(zt,nS,nT);
zts = (D1s*zt);
zts = zts(:);
zt = zt(:);
    
    
    
    
    
    
    
    
    
    
    
    
    
    