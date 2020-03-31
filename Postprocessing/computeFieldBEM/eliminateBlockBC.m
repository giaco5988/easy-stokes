%eliminate BC for block

function [fBCout,uBCout] = eliminateBlockBC(x,y,fBC,uBC,PARAM)

%compute location of the singularity
[~,~,nnn] = computeSingularityLocation(x,y,PARAM);

fBCout = fBC;
uBCout = uBC;
for i = 1:numel(PARAM.panels)
    
    %range
    [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,i);
    startMatrix = getSingRange(panelRange(1),nnn);
    [~,endMatrix] = getSingRange(panelRange(end),nnn);
    
   if PARAM.eliminateBlockPostprocessing(i)==1
       
       fBCout(2*startMatrix-1:2*endMatrix) = 0;
       uBCout(2*startMatrix-1:2*endMatrix) = 0;
       
   end
    
end