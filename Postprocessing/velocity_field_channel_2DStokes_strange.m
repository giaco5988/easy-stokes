%droplet VELOCITY FIELD visualization using uniform grid

%WON'T WORK BECAUSE BEM WILL TRY TO CONSERVE THE MASS WHILE THE BC REQUIRE A SOURCE!!!!!!!!!!!!!!!!!!! 
function [U,V,X,Y] = velocity_field_channel_2DStokes_strange(a,b,solution,resx,resy,PARAM)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

  n = PARAM.n;  m = PARAM.m;    j = PARAM.j;

  visc2 = PARAM.visc2;
  
  yy1 = max(b);
  yy2 = min(b);
  xx1 = min(a);
  xx2 = max(a);
  
  %mesh finess per unit length
  radial = round(resy*abs(yy1-yy2));
  axial = round(resx*abs(xx1-xx2));
    
    %velocities
    ux = [solution(1:2:2*n-1); zeros(m,1); solution(2*(n+m)+1:2:2*(n+m+j)-1); zeros(m,1)];
    uy = zeros(2*m+n+j,1);
    
    %stresses
    fx = [-PARAM.p_out*ones(n,1); solution(2*n+1:2:2*(n+m)-1,1); ...
        -(PARAM.press_grad+PARAM.p_out)*ones(j,1); solution(2*(n+m+j)+1:2:end-1)];
    fy = solution(2:2:end);

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(yy2,yy1,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu_2DStokes(a,b,X0,Y0);
    
    no = sqrt(diff(b).^2+diff(a).^2);
    r = [-diff(b)./no; diff(a)./no];
    
    R1 = repmat(r(1,1:end),numel(X0),1);
    R2 = repmat(r(2,1:end),numel(X0),1);
    
    A11 = TXXX.*R1+TXXY.*R2;
    A12 = TXYX.*R1+TXYY.*R2;
    A21 = TYXX.*R1+TYXY.*R2;
    A22 = TYYX.*R1+TYYY.*R2;
    
    %compute the integral to obtain the velocities
    U = (-GXX*fx-GXY*fy + visc2*(A11*ux+A12*uy))/pi/4;
    V = (-GYX*fx-GYY*fy + visc2*(A21*ux+A22*uy))/pi/4;
    
%     U = visc2*(A11*ux+A12*uy)/pi/4;
%     V = visc2*(A21*ux+A22*uy)/pi/4;
%     
%     U = (-GXX*fx-GXY*fy)/pi/4;
%     V = (-GYX*fx-GYY*fy)/pi/4;

    %reshape in the grid shape
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);

    figure
    quiver(X,Y,U,V)
    axis equal
    hold on
    quiver((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,ux',uy')
    if PARAM.continuity==1
        plot(PARAM.x_inj,PARAM.y_inj,'or')
    end
    hold off
    xlabel('x')
    ylabel('y')
    title(['flow field, \mu=' num2str(visc2)])
    
end