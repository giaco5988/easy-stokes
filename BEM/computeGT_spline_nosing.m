%compute Green's functions and associated stress tensor everywhere is
%needed with poz function assuming linear distribution of the variables
%here singular treatment is not necessary because the singularities are on
%the oyher bolck

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_nosing(xsing,ysing,ax,ay,bx,by,cx,cy,dx,dy)

    %tic

    %number of singularities
    N = numel(xsing);
    %number of element
    el = numel(ax);
    
    GXX = zeros(N);
    GXY = zeros(N);
    GYX = zeros(N);
    GYY = zeros(N);
    A11 = zeros(N);
    A12 = zeros(N);
    A21 = zeros(N);
    A22 = zeros(N);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax = reshape((repmat(ax,6,1)),1,6*el);
    bbx = reshape((repmat(bx,6,1)),1,6*el);
    ccx = reshape((repmat(cx,6,1)),1,6*el);
    ddx = reshape((repmat(dx,6,1)),1,6*el);
    aay = reshape((repmat(ay,6,1)),1,6*el);
    bby = reshape((repmat(by,6,1)),1,6*el);
    ccy = reshape((repmat(cy,6,1)),1,6*el);
    ddy = reshape((repmat(dy,6,1)),1,6*el);
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;
    
    phia = 1-GPt;
    phib = GPt;

    % point where the variable will be stored
    X0 = xsing;
    Y0 = ysing;
    
    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,el);
    
    %global coordianted retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %chabge sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    
%     figure
%     plot(eta,nx)
%     hold on
%     plot(eta,ny,'r')
%     hold off
    
%     figure
%     plot(eta,beta,'o-')
%     hold on
%     axis equal
%     plot(X0,Y0,'og-')
%     hold off
    
    h = sqrt(deta.*deta+dbeta.*dbeta);
    h0 = sqrt(bx.*bx+by.*by);
    h1 =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    manyh = repmat(h,N,1)';
    
%     figure
%     plot(x(1:end-1),h0,'o-')
%     hold on
%     plot(x(2:end),h1,'r-o')
%     plot(eta,h,'k-o')
%     hold off

%     figure
%     plot(h,'o-')

    X0matr = repmat(X0,6*el,1);
    Y0matr = repmat(Y0,6*el,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);

    [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
                sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
            
    SXX = manyh.*SXX;
    SXY = manyh.*SXY;
    SYX = manyh.*SYX;
    SYY = manyh.*SYY;
    QXXX = manyh.*QXXX;
    QXXY = manyh.*QXXY;
    QXYX = manyh.*QXYX;
    QXYY = manyh.*QXYY;
    QYXX = manyh.*QYXX;
    QYXY = manyh.*QYXY;
    QYYX = manyh.*QYYX;
    QYYY = manyh.*QYYY;
       
    %moltiplicate normals
    R1 = repmat(nx,N,1)';
    R2 = repmat(ny,N,1)';
    
    %for integration
    manyGW = repmat(GW',el,N);
    PHIA = repmat(phia',el,N);
    PHIB = repmat(phib',el,N);
    
    intSXXa = cumsum(SXX.*manyGW.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*PHIB);
    intSYXa = cumsum(SYX.*manyGW.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*PHIB);
    intQXXXa = cumsum(QXXX.*manyGW.*PHIA.*R1);
    intQXXXb = cumsum(QXXX.*manyGW.*PHIB.*R1);
    intQXXYa = cumsum(QXXY.*manyGW.*PHIA.*R2);
    intQXXYb = cumsum(QXXY.*manyGW.*PHIB.*R2);
    intQXYXa = cumsum(QXYX.*manyGW.*PHIA.*R1);
    intQXYXb = cumsum(QXYX.*manyGW.*PHIB.*R1);
    intQXYYa = cumsum(QXYY.*manyGW.*PHIA.*R2);
    intQXYYb = cumsum(QXYY.*manyGW.*PHIB.*R2);
    intQYXXa = cumsum(QYXX.*manyGW.*PHIA.*R1);
    intQYXXb = cumsum(QYXX.*manyGW.*PHIB.*R1);
    intQYXYa = cumsum(QYXY.*manyGW.*PHIA.*R2);
    intQYXYb = cumsum(QYXY.*manyGW.*PHIB.*R2);
    intQYYXa = cumsum(QYYX.*manyGW.*PHIA.*R1);
    intQYYXb = cumsum(QYYX.*manyGW.*PHIB.*R1);
    intQYYYa = cumsum(QYYY.*manyGW.*PHIA.*R2);
    intQYYYb = cumsum(QYYY.*manyGW.*PHIB.*R2);
    
    subSXXa(2:el,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:el,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:el,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:el,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:el,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:el,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:el,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:el,:) = intSYYb(6:6:end-6,:);
    subQXXXa(2:el,:) = intQXXXa(6:6:end-6,:);
    subQXXXb(2:el,:) = intQXXXb(6:6:end-6,:);
    subQXXYa(2:el,:) = intQXXYa(6:6:end-6,:);
    subQXXYb(2:el,:) = intQXXYb(6:6:end-6,:);
    subQXYXa(2:el,:) = intQXYXa(6:6:end-6,:);
    subQXYXb(2:el,:) = intQXYXb(6:6:end-6,:);
    subQXYYa(2:el,:) = intQXYYa(6:6:end-6,:);
    subQXYYb(2:el,:) = intQXYYb(6:6:end-6,:);
    subQYXXa(2:el,:) = intQYXXa(6:6:end-6,:);
    subQYXXb(2:el,:) = intQYXXb(6:6:end-6,:);
    subQYXYa(2:el,:) = intQYXYa(6:6:end-6,:);
    subQYXYb(2:el,:) = intQYXYb(6:6:end-6,:);
    subQYYXa(2:el,:) = intQYYXa(6:6:end-6,:);
    subQYYXb(2:el,:) = intQYYXb(6:6:end-6,:);
    subQYYYa(2:el,:) = intQYYYa(6:6:end-6,:);
    subQYYYb(2:el,:) = intQYYYb(6:6:end-6,:);
    
    %compute integral
    GXX(1:N,1:N-1) = (intSXXa(6:6:end,:)-subSXXa)';
    GXX(1:N,2:N) = GXX(1:N,2:N) + (intSXXb(6:6:end,:)-subSXXb)';
    GXY(1:N,1:N-1) = (intSXYa(6:6:end,:)-subSXYa)';
    GXY(1:N,2:N) = GXY(1:N,2:N) + (intSXYb(6:6:end,:)-subSXYb)';
    GYX(1:N,1:N-1) = (intSYXa(6:6:end,:)-subSYXa)';
    GYX(1:N,2:N) = GYX(1:N,2:N) + (intSYXb(6:6:end,:)-subSYXb)';
    GYY(1:N,1:N-1) = (intSYYa(6:6:end,:)-subSYYa)';
    GYY(1:N,2:N) = GYY(1:N,2:N) + (intSYYb(6:6:end,:)-subSYYb)';
    
    A11(:,1:N-1) = (intQXXXa(6:6:end,:)-subQXXXa + intQXXYa(6:6:end,:)-subQXXYa)';
    A11(:,2:N) = A11(:,2:N) + (intQXXXb(6:6:end,:)-subQXXXb + intQXXYb(6:6:end,:)-subQXXYb)';
    A12(:,1:N-1) = (intQXYXa(6:6:end,:)-subQXYXa + intQXYYa(6:6:end,:)-subQXYYa)';
    A12(:,2:N) = A12(:,2:N) + (intQXYXb(6:6:end,:)-subQXYXb + intQXYYb(6:6:end,:)-subQXYYb)';
    A21(:,1:N-1) = (intQYXXa(6:6:end,:)-subQYXXa + intQYXYa(6:6:end,:)-subQYXYa)';
    A21(:,2:N) = A21(:,2:N) + (intQYXXb(6:6:end,:)-subQYXXb + intQYXYb(6:6:end,:)-subQYXYb)';
    A22(:,1:N-1) = (intQYYXa(6:6:end,:)-subQYYXa + intQYYYa(6:6:end,:)-subQYYYa)';
    A22(:,2:N) = A22(:,2:N) + (intQYYXb(6:6:end,:)-subQYYXb + intQYYYb(6:6:end,:)-subQYYYb)';
    
    %T = toc

end