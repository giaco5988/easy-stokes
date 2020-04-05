%check metrics correlating physical soace and curvilinear coordinate

function [xNew,yNew] = makeUniformGrid(x,y,toleranceRemesh)

    TEMP.distr = [0];
    TEMP.normRemesh = [toleranceRemesh];
    TEMP.orderVariableStokes = [1];
    TEMP.adaptDistr = [];
    TEMP.orderGeometryStokes = [1];
    TEMP.SPlinesType = [2];
    TEMP.maxElem = 2e-2;
    TEMP.maxNumberTotalElem = 500;
    TEMP.minSizeElemRemesh = 1e-3;
    [xNew,yNew] = remeshDistrPanels2(1,x',y',TEMP,1);
    xNew = xNew';   yNew = yNew';
    yNew([1 end]) = [0 0];

end