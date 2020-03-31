%plot geoemtry

function plotGeometry(x,y,threeD,theta,color3D,transp,PARAM)

if threeD==1
    
    for i = 1:numel(PARAM.panels)

       PARAM.orderVariable = PARAM.orderVariableLaplace;
       PARAM.orderGeometry = PARAM.orderGeometryLaplace;

       %compute block coordinates
       [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);
       
       solidOfRevolution(yHere,xHere,theta{i},color3D{i},transp(i))
       hold on

    end
    
else

    for i = 1:numel(PARAM.panels)

       PARAM.orderVariable = PARAM.orderVariableLaplace;
       PARAM.orderGeometry = PARAM.orderGeometryLaplace;

       %compute block coordinates
       [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);

       if PARAM.fluxBC{i}==0
           col = 'k';
       else
           col = 'k';
       end

       if min(yHere)==0
            fill([xHere flip(xHere)],[yHere -flip(yHere)],'w')
       else
            fill(xHere,yHere,'w')
            fill(xHere,-yHere,'w')
       end

       plot(xHere,yHere,col)
       hold on
       plot(xHere,-yHere,col)

    end

end

axis equal