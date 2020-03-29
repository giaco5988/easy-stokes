%check metrics correlating physical soace and curvilinear coordinate

function activate = checkMetricCurvilinear(t,xyMode,PARAM)

activate = 0;
if (PARAM.remesh==1||PARAM.remesh==2) && t>=PARAM.firstTimeRemesh
    xMode = xyMode(1:2:end-1);
    yMode = xyMode(2:2:end);

    if PARAM.legendre==1||PARAM.legendre==2
        PPP = PARAM.PPP;
        x = LegendreBuildXY(xMode,PPP);
        y = LegendreBuildXY(yMode,PPP);
    elseif PARAM.legendre==0
        x = chebcoeffs2chebvals(xMode);
        y = chebcoeffs2chebvals(yMode);
    end

    %compute arclenght
    D1 = PARAM.D1;  w = PARAM.WG';
    l = computeTotalArcLengthSpectral(x,y,PARAM);
    xp = D1*x;    yp = D1*y;
    dl = sqrt(xp.^2+yp.^2)/l;
    
    if PARAM.normRemesh==1
        doRemesh = w*(dl-(D1*PARAM.remeshMapping(PARAM.t))).^2;
    elseif PARAM.normRemesh==2
        error('Not updated for mapping')
        ddl = D1*dl;
        doRemesh = w*ddl.^2;
    elseif PARAM.normRemesh==2
        doRemesh = w*(dl-(D1*SPECTRAL.remeshMapping(SPECTRAL.t))).^2 + w*dll.^2;
    end
    
    if doRemesh>PARAM.tolRemesh
    %if w*(dll/l).^2>PARAM.tolRemesh
        activate = 1;
        if PARAM.ODE==2||PARAM.ODE==9
            if PARAM.remesh==1 && PARAM.algorithm~=6
                display(['Remesh spectrally t=' num2str(t)])
            elseif PARAM.remesh==2 && PARAM.algorithm~=6
                display(['Remesh with splines t=' num2str(t)])
            end
        end
    end
    
end