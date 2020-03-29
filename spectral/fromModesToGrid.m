%compute modes coefficients from grid points in curvilinear coordinates

function [x,y] = fromModesToGrid(xMode,yMode,PARAM)

%compute current grid points
if PARAM.legendre==1||PARAM.legendre==2
        PPP = PARAM.PPP;
        x = LegendreBuildXY(xMode,PPP);
        y = LegendreBuildXY(yMode,PPP);
elseif PARAM.legendre==0
        addN = PARAM.n+1-numel(xMode);
        xMode = [xMode; zeros(addN,1)];
        yMode = [yMode; zeros(addN,1)];
        x = chebcoeffs2chebvals(xMode);
        y = chebcoeffs2chebvals(yMode);
end