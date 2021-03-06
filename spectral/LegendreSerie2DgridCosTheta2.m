%compute series of legendre form given data points

function f = LegendreSerie2DgridCosTheta2(theta,x,PPP,WG)
    
    %initialize variables
    [nS,nT] = size(x);
    modes = nT;
    manyWG = repmat(WG',nT,1);
    
    f = (PPP.*manyWG.*sin(theta))*x';
    coeff = repmat((0:modes-1)',1,nS)';
    f = f'.*(2*coeff+1)/2;

end