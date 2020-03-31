%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,X,Y] = velocity_field_channel_2DStokes_over(a,b,solution,resx,resy,PARAM)

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

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(yy2,yy1,radial);

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu_2DStokes(a,b,X0,Y0);
    %add velocity field and stresses with potential flow

    %in order to modify the physical boundary conditions
    [~,ux,uy,NXX,NXY] = gf_2D_Laplace((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,...
        PARAM.x_inj,PARAM.y_inj,PARAM.R,-PARAM.R,visc2,PARAM.p_out,PARAM.mass);
     
    %compute the velocities from potential flow that I will overlap
    [~,u,v] = gf_2D_Laplace(X0,Y0,...
        PARAM.x_inj,PARAM.y_inj,PARAM.R,-PARAM.R,visc2,PARAM.p_out,PARAM.mass);
    
    %normal vector
    no = sqrt(diff(b).^2+diff(a).^2);
    r = [-diff(b)./no; diff(a)./no];
    
    R1 = repmat(r(1,1:end),numel(X0),1);
    R2 = repmat(r(2,1:end),numel(X0),1);
    
    %double layer
    A11 = TXXX.*R1+TXXY.*R2;
    A12 = TXYX.*R1+TXYY.*R2;
    A21 = TYXX.*R1+TYXY.*R2;
    A22 = TYYX.*R1+TYYY.*R2;
    
    %stresses due to potential flow
    fx = NXX.*r(1,:)+NXY.*r(2,:);
    %fy = NYX.*r(1,:)+NYY.*r(2,:);
    
    %MODIFIED BOUNDARY CONDITIONS AND COMPLEMENTARY SOLUTION (from BEM)
    uux = [solution(1:2:2*n-1); -ux(n+1:m+n)'; solution(2*(n+m)+1:2:2*(n+m+j)-1); -ux(n+m+j+1:end)'];
    uuy = zeros(2*m+n+j,1);
    
    ffx = [-PARAM.p_out-fx(1:n)'; solution(2*n+1:2:2*(n+m)-1,1); ...
        -(PARAM.press_grad+PARAM.p_out)-fx(n+m+1:n+m+j)'; solution(2*(n+m+j)+1:2:end-1)];
    ffy = solution(2:2:end);
    
    %flow field from complmentary solution
    U = (-GXX*ffx-GXY*ffy + visc2*(A11*uux+A12*uuy))/pi/4;
    V = (-GYX*ffx-GYY*ffy + visc2*(A21*uux+A22*uuy))/pi/4;
    
    %physical flow field adding the potential
    U = U + u;
    V = V + v;

    %reshape in the grid shape
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);

    figure
    quiver(X,Y,U,V)
    axis equal
    hold on
    quiver((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,ux+uux',uy+uuy')
    if PARAM.continuity==1
        plot(PARAM.x_inj,PARAM.y_inj,'or')
    end
    hold off
    xlabel('x')
    ylabel('y')
    title(['flow field, \mu=' num2str(visc2)])
    
end