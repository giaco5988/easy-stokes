%compute droplet volume ic polar coordinate with spectral method

function V = VolumePolarAxisSpectral(theta,r,SPECTRAL)

    %derivatives and integration
    WG = SPECTRAL.WG;
    
    %compute volume
    V = 2/3*pi*WG'*(r.^3.*sin(theta));

end