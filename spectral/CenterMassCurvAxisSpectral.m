%compute droplet volume ic polar coordinate with spectral method

function xcm = CenterMassCurvAxisSpectral(x,y,SPECTRAL)

    %derivatives and integration
    WG = SPECTRAL.WG;
    D1 = SPECTRAL.D1;
    
    %derivatives
    xp = D1*x;
    
    %compute numerator
    num = WG'*(x.*y.^2.*xp);

    %compute denominator
    den = WG'*(y.^2.*xp);
    
    %compute center of mass
    xcm = num/den;
    
end