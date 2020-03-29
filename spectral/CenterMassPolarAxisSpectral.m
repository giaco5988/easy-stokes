%compute droplet volume ic polar coordinate with spectral method

function xcm = CenterMassPolarAxisSpectral(theta,r,SPECTRAL)

    %derivatives and integration
    WG = SPECTRAL.WG;
    
    %compute numerator
    num = 1/4*WG'*(r.^4.*cos(theta).*sin(theta));

    %compute denominator
    den = 1/3*WG'*(r.^3.*sin(theta));
    
    %compute center of mass
    xcm = num/den;
    
end