%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [PX,PY,PI1,PI2] = computeKernelPressureStokesAxisConst(x,y,x0,y0)

    %number of singularities
    N = numel(x0);
    elem = numel(x)-1;
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %compute normal vectors
    no = sqrt(diff(x).^2+diff(y).^2);
    nx = diff(y)./no;
    ny = -diff(x)./no;
    nnx = repmat(nx,6,1);
    nny = repmat(ny,6,1);
    nnx = repmat(nnx(:),1,N);
    nny = repmat(nny(:),1,N);
    
    %moltiplicate middle points coordinates
    xMiddle = (x(1:end-1)+x(2:end))/2;  yMiddle = (y(1:end-1)+y(2:end))/2;
    tempx = repmat(xMiddle,6,1);
    XX = reshape(tempx,1,6*elem);
    tempy = repmat(yMiddle,6,1);
    YY = reshape(tempy,1,6*elem);
    
    %lenght of the elements
    deltaX = diff(x);
    deltaY = diff(y);
    deltaL = sqrt(deltaX.^2+deltaY.^2);
    
    %every Gauss point
    GPX = repmat(GP,1,elem).*reshape((repmat(deltaX/2,6,1)),1,6*elem);
    GPY = repmat(GP,1,elem).*reshape((repmat(deltaY/2,6,1)),1,6*elem);
    
    globalX = XX+GPX;
    globalY = YY+GPY;
    
    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(globalX',1,N);
    globalYmatr = repmat(globalY',1,N);
        
    %compute gree function and the gradient
    [LX,LY,NXX,NXY,NYX,NYY] = kernelPressureAxis(globalXmatr,globalYmatr,X0matr,Y0matr);
      
    %multiply kernnel of doublelayer times the normal
    N1 = NXX.*nnx + NXY.*nny;
    N2 = NYX.*nnx + NYY.*nny;
    
    temp = reshape(repmat(deltaL/2,6,1),1,6*elem);
    manyDelta = repmat(temp',1,N);
    manyGW = repmat(GW',elem,N);
    
    %INTEGRATION
    intLX = cumsum(LX.*manyGW.*manyDelta);
    intLY = cumsum(LY.*manyGW.*manyDelta);
    intN1 = cumsum(N1.*manyGW.*manyDelta);
    intN2 = cumsum(N2.*manyGW.*manyDelta);
    
    subLX(2:elem,:) = intLX(6:6:end-6,:);
    subLY(2:elem,:) = intLY(6:6:end-6,:);
    subN1(2:elem,:) = intN1(6:6:end-6,:);
    subN2(2:elem,:) = intN2(6:6:end-6,:);
    
    PX = intLX(6:6:end,:)'-subLX';
    PY = intLY(6:6:end,:)'-subLY';
    PI1 = intN1(6:6:end,:)'-subN1';
    PI2 = intN2(6:6:end,:)'-subN2';

end