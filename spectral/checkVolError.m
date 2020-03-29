%check metrics correlating physical soace and curvilinear coordinate

function activate = checkVolError(xyMode,PARAM)

activate = 0;

%modes
xMode = xyMode(1:2:end-1);
yMode = xyMode(2:2:end);

%from modes to grid
[x,y] = fromModesToGrid(xMode,yMode,PARAM);

%compute volume
V = VolumeCurvilinearAxisSpectral(x,y,PARAM);
errV = abs((V-PARAM.V0))/PARAM.V0;

if errV>PARAM.volTol
    %disp('Execute volume correction')
    activate = 1;
end
    
end