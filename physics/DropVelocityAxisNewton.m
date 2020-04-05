%compute droplet velocity in axisymmetric domain

function R = DropVelocityAxisNewton(x,y,uNormal,Udrop)

    %center of mass
    xcm = center_mass(x,y);
    
    %function to integrate
    f = uNormal(Udrop).*(x'-xcm);
    
    %volume
    Volume = axis_int_gauss(x,y);
    
    %drop velocity as surface integral
    R = int_axis_spline_symmetric(x,y,f)/Volume;

end