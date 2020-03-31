%plot geoemtry

function plotGeometryStokes2d(x,y,PARAM)


    for i = 1:numel(PARAM.panels)

       PARAM.orderVariable = PARAM.orderVariableStokes;
       PARAM.orderGeometry = PARAM.orderGeometryStokes;

       %compute block coordinates
       [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);

       col = 'k-';

       plot(xHere,yHere,col)
       hold on
       plot(xHere,yHere,'.k','MarkerSize',18)

    end

axis equal