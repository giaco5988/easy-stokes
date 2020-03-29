%compute series of legendre form given data points

function f = LegendreSerieSpectralXY(x,PPP,SPECTRAL)

    %differentiation and integration
    WG = SPECTRAL.WG;
    
    %initialize variables
    modes = numel(x);

    f = sum(repmat(WG,1,modes).*repmat(x,1,modes).*PPP');
    coeff = 0:modes-1;
    f = f.*(2*coeff+1);
    f = f';
    

end