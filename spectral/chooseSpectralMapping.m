%choose spectral mapping

function mapping = chooseSpectralMapping(whichOne,intensity)

if whichOne==1    %normal mapping
    
    disp('No mapping, use spectral grid')
    
    mapping = @(s) s;
    
elseif whichOne==2
    
    disp(['Cluster point in the middle with intensity=' num2str(intensity)])
   
    intensity = intensity*pi;
    mapping = @(s) (tan((s-0.5)*intensity)/tan(intensity/2))/2+0.5;
    
elseif whichOne==3
    
    disp(['Cluster point to the right=' num2str(intensity)])
   
    %mapping = @(s) exp(s/(1-intensity))/exp(1/(1-intensity));
    mapping = @(s) s.^(1/(1-intensity));
    
    error('Not woking well, remesh fails')
    
else
    
    error('Not implemented')
    
end