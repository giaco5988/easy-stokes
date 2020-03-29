%compute droplet volume ic polar coordinate with spectral method

function A = surfaceCurvilinearAxisSpectral(x,y,SPECTRAL)

    %derivatives and integration
    WG = SPECTRAL.WG;
    D1 = SPECTRAL.D1;
    
    %compute volume
    xp = D1*x;  yp = D1*y;
    h = sqrt(xp.^2+yp.^2);
    A = 2*pi*WG'*(y.*h);

end