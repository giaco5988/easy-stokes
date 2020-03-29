%compute droplet velocity in axisymmetric domain

function Vdrop = DropVelocityAxisSpectralTheta(r,theta,uNormal,SPECTRAL)

    %integration weights
    WG = SPECTRAL.WG;

    %center of mass
    xcm = CenterMassPolarAxisSpectral(theta,r,SPECTRAL);
    
    %function to integrate
    f = 2*pi*uNormal.*(r.*cos(theta)-xcm).*r.*sin(theta);
    
    %metrics term
    rp = SPECTRAL.D1*r;
    h = sqrt(r.^2+rp.^2);
    
    %volume
    Volume = VolumePolarAxisSpectral(theta,r,SPECTRAL);
    
    %drop velocity as surface integral
    Vdrop = WG'*(f.*h)/Volume;

end