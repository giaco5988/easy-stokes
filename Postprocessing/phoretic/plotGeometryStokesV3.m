%plot geoemtry

function plotGeometryStokesV3(x,y,threeD,theta,color3D,transp,fillWhite,PARAM)

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
        
        if PARAM.eliminateBlockPostprocessing(i)==0

           PARAM.orderVariable = PARAM.orderVariableStokes;
           PARAM.orderGeometry = PARAM.orderGeometryStokes;

           %compute block coordinates
           [xHere,yHere] = getBlockCoordinates(x,y,PARAM,i);

           col = 'k';

           if min(yHere)<1e-15 && fillWhite(i)==1
                fill([xHere flip(xHere)],[yHere -flip(yHere)],'w')
           elseif min(yHere)<1e-15 && fillWhite(i)==2
                fill([xHere flip(xHere)],[yHere -flip(yHere)],'k')
           elseif fillWhite(i)==1
                fill(xHere,yHere,'w')
                fill(xHere,-yHere,'w')
           elseif fillWhite(i)==2
                fill(xHere,yHere,'k')
                fill(xHere,-yHere,'k')
           end

           plot(xHere,yHere,col)
           hold on
           plot(xHere,-yHere,col)
       
        end

    end

end

%axis equal