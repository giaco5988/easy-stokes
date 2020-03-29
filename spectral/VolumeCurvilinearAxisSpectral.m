%compute droplet volume ic polar coordinate with spectral method

function V = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL)

    %derivatives and integration
    WG = SPECTRAL.WG;
    D1 = SPECTRAL.D1;
    
    %compute volume
    xp = D1*x;
    V = -pi*WG'*(y.^2.*xp);

end