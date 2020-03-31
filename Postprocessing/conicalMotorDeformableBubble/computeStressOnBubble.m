%stress on bubble for Post-Processing

function [fBCx,fBCy,dfXrep,dfYrep] = computeStressOnBubble(x,y,panel,PARAM)

%compute curvature
[nx,ny] = computeNormalVector(x{panel},y{panel},PARAM.orderVariable(panel),PARAM.orderGeometry(panel),PARAM.SPlinesType(panel));
[k1,k2] = computeCurvatureSplines(x{panel},y{panel},PARAM.orderVariable(panel));
k = k1+k2;
          
%known BC
fBCx = PARAM.stressBC{panel}*k.*nx;
fBCy = PARAM.stressBC{panel}*k.*ny;

%add repulsive force for droplet (or bubble)
for i = 1:numel(PARAM.panels)
      
      if (PARAM.repulsiveForces(i)==3||PARAM.repulsiveForces(i)==4) && PARAM.blockType(i)==2
      
      if numel(PARAM.panels)>2
          error('Current implementation works only for two block')
      end
      
      for k = 1:numel(PARAM.panels)
          
            if k~=i
      
                [dfXrep,dfYrep] = disjoiningPressureBlocks(x,y,i,k,PARAM);
            
            end
            
      end
    
      end
    
end