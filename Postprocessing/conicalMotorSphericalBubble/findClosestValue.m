

function valueOut = findClosestValue(xLine,yLine,X,Y,valueIN)

valueOut = zeros(numel(xLine),1);
for i = 1:numel(xLine)
    
    [~,indX] = min(abs(xLine(i)-X));
    [~,indY] = min(abs(yLine(i)-Y));
    valueOut(i) = valueIN(indY,indX);
    
end