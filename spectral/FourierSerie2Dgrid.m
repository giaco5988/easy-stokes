%compute series of legendre form given data points

function [an,bn] = FourierSerie2Dgrid(x,WG,s)
    
    %initialize variables
    [nS,~] = size(x);
    manyWG = repmat(WG',nS,1);
    
    n = 0:nS-1;
    nnn = repmat(n',1,nS);
    manyCos = cos(repmat(s',nS,1).*nnn);
    manySin = sin(repmat(s',nS,1).*nnn);
    
    an = (manyCos.*manyWG)*x/pi;
    bn = (manySin.*manyWG)*x/pi;
    
    %dealiasing
    N = nS/2;
    an(N+1:end,:) = 0;
    bn(N+1:end,:) = 0;

end