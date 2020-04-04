%compute function which gives the modification of the volume due to the
%motion of one point

function R = ModifyVolumeXYmodes(xBase,yBase,nx,ny,DN,V0,perturbMode)

    perturbMode = [DN; perturbMode];
    perturb = chebcoeffs2chebvals(perturbMode);
    perturbCheb = chebfun(perturb,[0 1]);
    perturb = perturbCheb(linspace(0,1,numel(xBase)));
    
    %displace in the normal direction
    x = xBase + perturb'.*nx';
    y = yBase + perturb'.*ny';
    
    %residual
    R = axis_int_gauss_vect(x',y')-V0;

end