%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeXY(x,y,nx,ny,DN,V0)
    
    %displace in the normal direction
    x = x + DN*nx;
    y = y + DN*ny;
    
    %residual
    R = axis_int_gauss_vect(x,y)-V0;

end