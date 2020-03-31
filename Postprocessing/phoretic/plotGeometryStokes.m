%plot geoemtry

function plotGeometryStokes(x,y,threeD,theta,color3D,transp,fillWhite,PARAM)

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
            fill([xHere flip(xHere)],[yHere -flip(yHere)],'w')
       elseif fillWhite==1
            fill(xHere,yHere,'w')
            fill(xHere,-yHere,'w')
       end

       plot(xHere,yHere,col)
       hold on
       plot(xHere,-yHere,col)

    end

end

axis equal