%plot geoemtry

function plotGeometryDropRemesh(x,y,PARAM)

for i = 1:numel(PARAM.n)
    
   if PARAM.velBC{i}==0
       col = 'kx-';
   else
       col = 'xr-';
   end
    
   xHere = x{i};
   yHere = y{i};
   
   plot(xHere,yHere,col)
   hold on
   plot(xHere,-yHere,col)
    
end