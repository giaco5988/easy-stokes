%plot geoemtry

function plotGeometryStokesV2(x,y,threeD,theta,color3D,transp,fillWhite,shiftUP,PARAM)

if threeD==1
    
    for i = 1:numel(PARAM.panels)

       PARAM.orderVariable = PARAM.orderVariableStokes;
       PARAM.orderGeometry = PARAM.orderGeometryStokes;

       %compute block coordinates
       [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);
       
       solidOfRevolution(yHere,xHere,theta{i},color3D{i},transp(i))
       hold on

    end

    
    
else

    for i = 1:numel(PARAM.panels)

       PARAM.orderVariable = PARAM.orderVariableStokes;
       PARAM.orderGeometry = PARAM.orderGeometryStokes;

       %compute block coordinates
       [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);

       col = 'k';

       if min(yHere)<1e-15 && fillWhite
            fill([xHere flip(xHere)],[yHere -flip(yHere)]+shiftUP,'w')
       elseif fillWhite==1
            fill(xHere,yHere+shiftUP,'w')
            fill(xHere,-yHere+shiftUP,'w')
       end

       plot(xHere,yHere+shiftUP,col)
       hold on
       plot(xHere,-yHere+shiftUP,col)

    end

end

axis equal