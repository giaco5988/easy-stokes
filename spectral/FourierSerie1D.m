%compute series of legendre form given data points

function [an,bn] = FourierSerie1D(x,f,WG)
    
    mult = 1/pi;
    
    %initialize
    an = zeros(numel(f),1);
    bn = zeros(numel(f),1);
    
    %compute coefficients
    an(1) = WG'*f*mult;
    for i = 2:numel(f)
        
        n = i-1;
        an(i) = WG'*(f.*cos(n*x))*mult;
        bn(i) = WG'*(f.*sin(n*x))*mult;
        
    end
    
    N = numel(x)/2;
    an(N+1:end) = 0;
    bn(N+1:end) = 0;

end