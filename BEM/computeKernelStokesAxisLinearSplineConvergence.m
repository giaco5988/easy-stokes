%compute Green's functions and associated stress tensor everywhere is
%needed with poz function

function [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeKernelStokesAxisLinearSplineConvergence(x0,y0,indIn,indOut,ST,ax,ay,bx,by,cx,cy,dx,dy,kernelFreeSpace,wallPos)

    %number of singularities
    N = numel(x0);
    elem = numel(ax);
    
    %initialize
    GXX = zeros(N,elem+1);
    GXY = zeros(N,elem+1);
    GYX = zeros(N,elem+1);
    GYY = zeros(N,elem+1);
    A11 = zeros(N,elem+1);
    A12 = zeros(N,elem+1);
    A21 = zeros(N,elem+1);
    A22 = zeros(N,elem+1);
    
    %gauss points an weigths
    GP = [-0.932469514203152 -0.661209386466265 -0.238619186083197 0.238619186083197 0.661209386466265 0.932469514203152];
    GW = [0.171324492379170 0.360761573048139 0.467913934572691 0.467913934572691 0.360761573048139 0.171324492379170];
    
    aax = reshape((repmat(ax,6,1)),1,6*elem);
    bbx = reshape((repmat(bx,6,1)),1,6*elem);
    ccx = reshape((repmat(cx,6,1)),1,6*elem);
    ddx = reshape((repmat(dx,6,1)),1,6*elem);
    aay = reshape((repmat(ay,6,1)),1,6*elem);
    bby = reshape((repmat(by,6,1)),1,6*elem);
    ccy = reshape((repmat(cy,6,1)),1,6*elem);
    ddy = reshape((repmat(dy,6,1)),1,6*elem);    
    
    % transform GP because the spline parammeter is defined between 0 and 1
    GPt = (GP+1)/2;
    GW = GW/2;
    
    phia = 1-GPt;
    phib = GPt;

    % point where the variable will be stored
    X0 = x0';
    Y0 = y0';
    
    %points where I perform gauss integration
    GPglobal = repmat(GPt,1,elem);
    
    %global coordianted retrived on the splines
    eta = aax+bbx.*GPglobal+ccx.*GPglobal.*GPglobal+ddx.*GPglobal.*GPglobal.*GPglobal;
    beta = aay+bby.*GPglobal+ccy.*GPglobal.*GPglobal+ddy.*GPglobal.*GPglobal.*GPglobal;
    deta = bbx+2*ccx.*GPglobal+3*ddx.*GPglobal.*GPglobal;
    dbeta = bby+2*ccy.*GPglobal+3*ddy.*GPglobal.*GPglobal;
    
    %chabge sign for my convention
    nx = dbeta./sqrt(deta.*deta+dbeta.*dbeta);
    ny = -deta./sqrt(deta.*deta+dbeta.*dbeta);
    xp0 = bx;
    yp0 = by;
    xp1 = bx+2*cx+3*dx;
    yp1 = by+2*cy+3*dy;
    
    %moltiplicate normals
    nnx = repmat(nx,N,1)';
    nny = repmat(ny,N,1)';
    
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

%     h = zeros(1,el*6);
%     h0 = zeros(el,1);
%     h1 = zeros(el,1);
    
    hTemp = sqrt(deta.*deta+dbeta.*dbeta);
    hx = hTemp.*nx;
    hy = hTemp.*ny;
    h0temp = sqrt(bx.*bx+by.*by);
    nx0 = yp0./h0temp;
    h0x = h0temp.*nx0;
    ny0 = -xp0./h0temp;
    h0y = h0temp.*ny0;
    h1temp =  sqrt((bx+2*cx+3*dx).*(bx+2*cx+3*dx)+(by+2*cy+3*dy).*(by+2*cy+3*dy));
    nx1 = yp1./h1temp;
    h1x = h1temp.*nx1;
    ny1 = -xp1./h1temp;
    h1y = h1temp.*ny1;
    manyh = repmat(hTemp,N,1)';
    
%     figure
%     plot(x(1:end-1),h0,'o-')
%     hold on
%     plot(x(2:end),h1,'r-o')
%     plot(eta,h,'k-o')
%     hold off

%     figure
%     plot(h,'o-')

    X0matr = repmat(X0,6*elem,1);
    Y0matr = repmat(Y0,6*elem,1);
    
    globalXmatr = repmat(eta',1,N);
    globalYmatr = repmat(beta',1,N);
        
    %compute gree function and the gradient
    if kernelFreeSpace==1     %free space
          [SXX,SXY,SYX,SYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] =...
            sgf_ax_fs_vect3 (globalXmatr,globalYmatr,X0matr,Y0matr);
    elseif kernelFreeSpace==2
            
            [SXX,SXY,SYX,SYY ...
            ,QXXX,QXXY,QXYX,QXYY ...
            ,QYXX,QYXY,QYYX,QYYY] = sgf_ax_w_vect (globalXmatr,globalYmatr,X0matr,Y0matr,wallPos);
        
    else
        error('Not implemented')
    end
    
    %mutuply kernnel of doublelayer times the normal
    manyhX = manyh.*nnx;
    manyhY = manyh.*nny;
    SXX = manyhX.*SXX;
    SXY = manyhY.*SXY;
    SYX = manyhX.*SYX;
    SYY = manyhY.*SYY;
    Q11 = manyhX.*QXXX + manyhY.*QXXY;
    Q12 = manyhX.*QXYX + manyhY.*QXYY;
    Q21 = manyhX.*QYXX + manyhY.*QYXY;
    Q22 = manyhX.*QYYX + manyhY.*QYYY;
      
    %singular treatment
    if ST==1
        
        range = indIn:indOut;
        [SXX(:,range),SYY(:,range),anA1x,anB1x,anA2x,anB2x,anA1y,anB1y,anA2y,anB2y] = stokesSTlinearSplineConvergence(SXX(:,range),SYY(:,range),eta,beta,hx,hy,h0x,h0y,h1x,h1y,x0(range),y0(range),GPt);
          
        
    end
    
    %INTEGRATION
    manyGW = repmat(GW',elem,N);
    PHIA = repmat(phia',elem,N);
    PHIB = repmat(phib',elem,N);
    
    %INTEGRATION
    intSXXa = cumsum(SXX.*manyGW.*PHIA);
    intSXXb = cumsum(SXX.*manyGW.*PHIB);
    intSXYa = cumsum(SXY.*manyGW.*PHIA);
    intSXYb = cumsum(SXY.*manyGW.*PHIB);
    intSYXa = cumsum(SYX.*manyGW.*PHIA);
    intSYXb = cumsum(SYX.*manyGW.*PHIB);
    intSYYa = cumsum(SYY.*manyGW.*PHIA);
    intSYYb = cumsum(SYY.*manyGW.*PHIB);
    intQ11a = cumsum(Q11.*manyGW.*PHIA);
    intQ11b = cumsum(Q11.*manyGW.*PHIB);
    intQ12a = cumsum(Q12.*manyGW.*PHIA);
    intQ12b = cumsum(Q12.*manyGW.*PHIB);
    intQ21a = cumsum(Q21.*manyGW.*PHIA);
    intQ21b = cumsum(Q21.*manyGW.*PHIB);
    intQ22a = cumsum(Q22.*manyGW.*PHIA);
    intQ22b = cumsum(Q22.*manyGW.*PHIB);
    
    subSXXa(2:elem,:) = intSXXa(6:6:end-6,:);
    subSXXb(2:elem,:) = intSXXb(6:6:end-6,:);
    subSXYa(2:elem,:) = intSXYa(6:6:end-6,:);
    subSXYb(2:elem,:) = intSXYb(6:6:end-6,:);
    subSYXa(2:elem,:) = intSYXa(6:6:end-6,:);
    subSYXb(2:elem,:) = intSYXb(6:6:end-6,:);
    subSYYa(2:elem,:) = intSYYa(6:6:end-6,:);
    subSYYb(2:elem,:) = intSYYb(6:6:end-6,:);
    subQ11a(2:elem,:) = intQ11a(6:6:end-6,:);
    subQ11b(2:elem,:) = intQ11b(6:6:end-6,:);
    subQ12a(2:elem,:) = intQ12a(6:6:end-6,:);
    subQ12b(2:elem,:) = intQ12b(6:6:end-6,:);
    subQ21a(2:elem,:) = intQ21a(6:6:end-6,:);
    subQ21b(2:elem,:) = intQ21b(6:6:end-6,:);
    subQ22a(2:elem,:) = intQ22a(6:6:end-6,:);
    subQ22b(2:elem,:) = intQ22b(6:6:end-6,:);
    
    %compute integral
    GXX(:,1:elem) = (intSXXa(6:6:end,:)-subSXXa)'; %+ vertcat(zeros(1,N-1),diag(-deltaL(1:N-1).*log(deltaL(1:N-1))+deltaL(1:N-1)));
    GXX(:,2:elem+1) = GXX(:,2:elem+1) + intSXXb(6:6:end,:)'-subSXXb';
    GXY(:,1:elem) = intSXYa(6:6:end,:)'-subSXYa';
    GXY(:,2:elem+1) = GXY(:,2:elem+1) + intSXYb(6:6:end,:)'-subSXYb';
    GYX(:,1:elem) = intSYXa(6:6:end,:)'-subSYXa';
    GYX(:,2:elem+1) = GYX(:,2:elem+1) + intSYXb(6:6:end,:)'-subSYXb';
    GYY(:,1:elem) = intSYYa(6:6:end,:)'-subSYYa';
    GYY(:,2:elem+1) = GYY(:,2:elem+1) + intSYYb(6:6:end,:)'-subSYYb';
    
    A11(:,1:elem) = (intQ11a(6:6:end,:)-subQ11a)';
    A11(:,2:elem+1) = A11(:,2:elem+1) + (intQ11b(6:6:end,:)-subQ11b)';
    A12(:,1:elem) = (intQ12a(6:6:end,:)-subQ12a)';
    A12(:,2:elem+1) = A12(:,2:elem+1) + (intQ12b(6:6:end,:)-subQ12b)';
    A21(:,1:elem) = (intQ21a(6:6:end,:)-subQ21a)';
    A21(:,2:elem+1) = A21(:,2:elem+1) + (intQ21b(6:6:end,:)-subQ21b)';
    A22(:,1:elem) = (intQ22a(6:6:end,:)-subQ22a)';
    A22(:,2:elem+1) = A22(:,2:elem+1) + (intQ22b(6:6:end,:)-subQ22b)';
    
    %add anaytical integration of the singularity
    if ST==1
        
        N = numel(h0x)+1;
        
        range1 = (2:N-2)+indIn(1)-1;
        range2 = (2:N-2)+indIn(1);
        
        %add analtyical integration
        GXX(range1,2:N-2) = GXX(range1,2:N-2) + anA1x;
        GXX(range1,3:N-1) = GXX(range1,3:N-1) + anB1x;
        GXX(range2,2:N-2) = GXX(range2,2:N-2) + anB2x;
        GXX(range2,3:N-1) = GXX(range2,3:N-1) + anA2x;
        
        GYY(range1,2:N-2) = GYY(range1,2:N-2) + anA1y;
        GYY(range1,3:N-1) = GYY(range1,3:N-1) + anB1y;
        GYY(range2,2:N-2) = GYY(range2,2:N-2) + anB2y;
        GYY(range2,3:N-1) = GYY(range2,3:N-1) + anA2y;
        
        %because singularity on the axis don't has to be treated
        GXX(range(2),1) = GXX(range(2),1) + 1/2*h1x(1);
        GXX(range(2),2) = GXX(range(2),2) + 3/2*h1x(1);
        GXX(range(end-1),end-1) = GXX(range(end-1),end-1) + 3/2*h0x(elem);
        GXX(range(end-1),end) = GXX(range(end-1),end) + 1/2*h0x(elem);
        
        GYY(range(2),1) = GYY(range(2),1) + 1/2*h1y(1);
        GYY(range(2),2) = GYY(range(2),2) + 3/2*h1y(1);
        GYY(range(end-1),end-1) = GYY(range(end-1),end-1) + 3/2*h0y(elem);
        GYY(range(end-1),end) = GYY(range(end-1),end) + 1/2*h0y(elem);
    
    end

end