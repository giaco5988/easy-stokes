

function [valueOut,outFit] = findClosestValueInterp(xLine,yLine,X,Y,valueIN)

[Xgrid,Ygrid] = ndgrid(X,Y);
valueIN = valueIN';
outFit = fit([Xgrid(:), Ygrid(:)],valueIN(:),'linearinterp');
valueOut = outFit(xLine,yLine);