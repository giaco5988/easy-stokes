%add repuslive forces in Post Processing

function fBC = addRepulsivePostprocessing(x,y,nnn,fBC,PARAM)

for i = 1:numel(PARAM.panels)
      
      if (PARAM.repulsiveForces(i)==3||PARAM.repulsiveForces(i)==4) && PARAM.blockType(i)==2
      
      if numel(PARAM.panels)>2
          error('Current implementation works only for two block')
      end
      
      %range
      [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,i);
      startMatrix = getSingRange(panelRange(1),nnn);
      [~,endMatrix] = getSingRange(panelRange(end),nnn);
      
      for k = 1:numel(PARAM.panels)
          
            if k~=i
      
                [dfXrep,dfYrep] = disjoiningPressureBlocks(x,y,i,k,PARAM);
                fBC(2*startMatrix:2:2*endMatrix) = fBC(2*startMatrix:2:2*endMatrix) + dfYrep;
                if PARAM.repulsiveForces(i)==4
                   fBC(2*startMatrix-1:2:2*endMatrix-1) = fBC(2*startMatrix-1:2:2*endMatrix-1) + dfXrep;
                end
            
            end
            
      end
    
      end
    
end