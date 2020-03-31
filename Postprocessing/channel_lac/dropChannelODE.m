%post processinf drop in pipe from ODE

close all

%options
plotShape = 1;
step = 10;

%number of iteration
ite = round(numel(T)/step);

%initialize
V = zeros(1,ite);
D = zeros(1,ite);
V0 = 4/3*pi*PARAM.alpha^3;

for k = 1:ite
    
    i = 1+(k-1)*step;
    
    display(['loop' num2str(k) ' of ' num2str(ite)])
    
    %coordinates at this operation
    xDrop = Y(i,1:2:end-1);     yDrop = Y(i,2:2:end);
    
    %compute location of ghost point
    fVolume = @(rho) ModifyVolume4(xDrop,yDrop,rho,V0);
    options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
    xEnd = fsolve(fVolume,xDrop(end),options);
    
    %compute drop coordinates with ghost point
    xDrop = [xDrop xEnd];   yDrop = [yDrop 0];
    
    %compute droplet volume
    V(k) = axis_int_gauss_vect(xDrop,yDrop);
    
    %check shape
    aaa = max(xDrop)-min(xDrop);
    bbb = max(yDrop);
    D(k) = (aaa-bbb)/(aaa+bbb);
    
    %plot shape
    if plotShape==1
       
       subplot(3,1,1)
       plot(xDrop,yDrop,'r')
       hold on
       plot(xDrop,-yDrop,'r')
       plot(xWall,yWall,'k')
       plot(xWall,-yWall,'k')
       axis equal
       grid on
       axis([min(xDrop)-1 max(xDrop)+1 -max(yWall)-0.1 max(yWall)+0.1])
       xlabel('x')
       ylabel('r')
       hold off
       title('Droplet in channel')
       drawnow
        
    end
    
end

subplot(3,1,2)
errV = (V-V0)/V0;
plot(T(1:step:ite*step),errV,'o-')
grid on
xlabel('t')
ylabel('errV')
title('error on volume')

subplot(3,1,3)
plot(T(1:step:ite*step),D,'o-')
grid on
xlabel('t')
ylabel('D')
title('Deformation')












