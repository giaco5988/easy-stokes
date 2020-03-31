%when computing quantities in the domain with BEM, move point really close
%to the interface on the closest node of the interface and assign the
%corrisponding quantity of the node

function [Xsing,Ysing,ux,uy,p,InOut] = closeInterfaceVisuBlocks(x,y,xStore,yStore,Xsing,Ysing,ux,uy,p,InOut,uBC,maxElem,coeffSub,velBlocks,PARAM)

    %compute number of sing per panel
    [~,~,nnn] = computeSingularityLocation(x,y,PARAM);

    %transform normal-tangent basis in (z,r) basis
    for k = 1:numel(PARAM.n)
        
       if PARAM.typeBCstokes(k)==5
           
          [nx,ny] = computeNormalVector(x{k},y{k},PARAM.orderVariable(k),PARAM.orderGeometry(k),PARAM.SPlinesType(k)); 
          [startMatrix,endMatrix] = getSingRange(k,nnn);
          uBC(2*startMatrix-1:2:2*endMatrix-1) = uBC(2*startMatrix-1:2:2*endMatrix-1).*nx';
          uBC(2*startMatrix:2:2*endMatrix) = uBC(2*startMatrix:2:2*endMatrix).*ny';
          
       end
        
    end

    %add velocity BC for moving walls
    for k = 1:numel(PARAM.panels)
       
        %range
        [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,k);
        startMatrix = getSingRange(panelRange(1),nnn);
        [~,endMatrix] = getSingRange(panelRange(end),nnn);
        
        if PARAM.blockType(k)==1
            
            uBC(2*startMatrix-1:2:2*endMatrix-1) = uBC(2*startMatrix-1:2:2*endMatrix-1)+velBlocks(k);
            
        end
        
    end
    
    [Xsing,Ysing,ux,uy,p,InOut] = closeInterfaceVisu(xStore,yStore,Xsing,Ysing,ux,uy,p,InOut,uBC,maxElem,coeffSub);
    
end