%compute series of legendre form given data points

function f = LegendreSerie2Dgrid(x,PPP,WG)
    
    %initialize variables
    [nS,nT] = size(x);
    modes = nT;
    manyWG = repmat(WG',nT,1);
    
    f = (PPP.*manyWG)*x';
    coeff = repmat((0:modes-1)',1,nS)';
    f = f'.*(2*coeff+1)/pi;

end