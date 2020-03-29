%compute droplet velocity in axisymmetric domain

function Vdrop = DropVelocityAxisCurvilinear(x,y,uNormal,SPECTRAL)

    %center of mass
    xcm = CenterMassCurvAxisSpectral(x,y,SPECTRAL);
    
    %function to integrate
    f = uNormal.*(x-xcm);
    
    %volume
    Volume = VolumeCurvilinearAxisSpectral(x,y,SPECTRAL);
    
    %drop velocity as surface integral
    Vdrop = axisIntSpectralCurvilinear(x,y,f,SPECTRAL)/Volume;

end