%plot geoemtry

function plotGeometryDrop(x,y,PARAM,plotMesh)

for i = 1:numel(PARAM.n)
    
   if PARAM.velBC{i}==0
       col = 'k';
       colMesh = 'xk';
   else
       col = 'r';
       colMesh = 'rx';
   end
    
   xHere = x{i};
   yHere = y{i};
   
   plot(xHere,yHere,col)
   hold on
   if plotMesh==1
       plot(xHere,-yHere,colMesh)
   else
       plot(xHere,-yHere,col)
   end
   axis equal
    
end