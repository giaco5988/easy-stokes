%dealiasing point on the physical grid

function [xNew,yNew] = dealiasingGridXY(x,y,PARAM)

%compute the modes from the grid points
if PARAM.legendre==1||PARAM.legendre==2
       
       PPP = PARAM.PPP;
       xMode = LegendreSerieSpectralXY(x,PPP,PARAM);
       yMode = LegendreSerieSpectralXY(y,PPP,PARAM);
       
elseif PARAM.legendre==0
       
       x = chebfun(x,[0 1]);     y = chebfun(y,[0 1]);
       xMode = chebcoeffs(x);
       yMode = chebcoeffs(y);
       
       if numel(xMode)<PARAM.dealiasing
           xMode = [xMode; zeros(PARAM.dealiasing-numel(xMode),1)];
           yMode = [yMode; zeros(PARAM.dealiasing-numel(yMode),1)];
       end
   
end

%perform dealiasing
xMode(PARAM.dealiasing+1:end) = 0;
yMode(PARAM.dealiasing+1:end) = 0;

%compute de-aliased grid points
if PARAM.legendre==1||PARAM.legendre==2
        
        xNew = LegendreBuildXY(xMode,PPP);
        yNew = LegendreBuildXY(yMode,PPP);
        
elseif PARAM.legendre==0
        
        xNew = chebcoeffs2chebvals(xMode);
        yNew = chebcoeffs2chebvals(yMode);
        
end