%choose non linear function, 1 drop, Spectral

function f = tutorial_functionDropSpectral(PARAM)

if PARAM.legendre==1||PARAM.legendre==2
        f = @(t,xyMode) tutorial_dropLegendreCurvilinearModes(t,xyMode,PARAM);
elseif PARAM.legendre==0
        f = @(t,xyMode) tutorial_dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
end