%droplet VELOCITY FIELD visualization using uniform grid

function [U,V,p,X,Y] = velocity_2DStokes_inner(a,b,solution,resx,resy,PARAM)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

  n = PARAM.n;  m = PARAM.m;    j = PARAM.j;
  
  fixed = n+m+j+m;

  visc2 = PARAM.visc2;
  
  yy1 = max(b)-0.1;
  yy2 = min(b)+0.1;
  xx1 = min(a)+0.1;
  xx2 = max(a)-0.1;
  
  yy1 = max(b)-0.2;
  yy2 = min(b)+0.2;
  xx1 = min(a)+4;
  xx2 = max(a)-4;
  
  %mesh finess per unit length
  radial = round(resy*abs(yy1-yy2));
  axial = round(resx*abs(xx1-xx2));

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(yy2,yy1,radial);
    
    Xsing = [(a(1:fixed)+a(2:fixed+1))/2 (a(fixed+2:end-1)+a(fixed+3:end))/2];
    Ysing = [(b(1:fixed)+b(2:fixed+1))/2 (b(fixed+2:end-1)+b(fixed+3:end))/2];

    [X,Y] = meshgrid(x,y);
    
    X0 = X(:);
    Y0 = Y(:);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu_2DStokes_interface(a,b,X0,Y0,fixed);
    %add velocity field and stresses with potential flow
    
    %compute GF for pressure
    [PX,PY,PIXX,PIXY,PIYX,PIYY] = computeP_visu_2DStokes_interface(a,b,X0,Y0,fixed);
    
    %normal vector
    no = [sqrt(diff(b(1:fixed+1)).^2+diff(a(1:fixed+1)).^2) ...
      sqrt(diff(b(fixed+2:end)).^2+diff(a(fixed+2:end)).^2)];
    r = [-diff(b(1:fixed+1))./no(1:fixed) diff(b(fixed+2:end))./no(fixed+1:end);...
        diff(a(1:fixed+1))./no(1:fixed) -diff(a(fixed+2:end))./no(fixed+1:end)];
    
    t = [-diff(a(fixed+2:end))./no(fixed+1:end); -diff(b(fixed+2:end))./no(fixed+1:end)];
    
    %velocities from solution and BC
    u = PARAM.mass/2/pi/PARAM.alpha/PARAM.R;  %velocity in every point
    ux = u*r(1,n+m+m+j+1:end);
    uy = u*r(2,n+m+m+j+1:end);
    %ux = u*t(1,:);
    %uy = u*t(2,:);
    if PARAM.stress_in==1
        ux = [solution(1:2:2*n-1); zeros(m,1); solution(2*(n+m)+1:2:2*(n+m+j)-1); zeros(m,1); ux'];
        uy = [zeros(n+m+m+j,1); uy'];
    elseif PARAM.stress_in==0
        ux = [solution(1:2:2*n-1); zeros(m,1); 3/4*PARAM.Q/PARAM.R*(1-(Ysing(n+m+1:n+m+j)'/PARAM.R).^2); zeros(m,1); ux'];
        uy = [zeros(n+m+m+j,1); uy'];
    end
        
    
    if PARAM.stress_in==1
        %stresses from solution and BC
        fx = [PARAM.p_out*ones(n,1); solution(2*n+1:2:2*(n+m)-1,1); ...
            -(PARAM.press_grad+PARAM.p_out)*(ones(j,1)); solution(2*(n+m+j)+1:2:end-1)];
        fy = solution(2:2:end);
    elseif PARAM.stress_in==0
        %stresses from solution and BC
        fx = [PARAM.p_out*ones(n,1); solution(2*n+1:2:2*(n+m)-1,1); ...
            solution(2*(n+m)+1:2:2*(n+m+j)-1); solution(2*(n+m+j)+1:2:end-1)];
        fy = solution(2:2:end);
    end
        
    
    %rep normals
    R1 = repmat(r(1,1:end),numel(X0),1);
    R2 = repmat(r(2,1:end),numel(X0),1);
    
    %double layer
    A11 = TXXX.*R1+TXXY.*R2;
    A12 = TXYX.*R1+TXYY.*R2;
    A21 = TYXX.*R1+TYXY.*R2;
    A22 = TYYX.*R1+TYYY.*R2;
    
    %double layer pressure
    B1 = PIXX.*R1+PIXY.*R2;
    B2 = PIYX.*R1+PIYY.*R2;
    
    %flow field
    U = (-GXX*fx-GXY*fy + visc2*(A11*ux+A12*uy))/pi/4;
    V = (-GYX*fx-GYY*fy + visc2*(A21*ux+A22*uy))/pi/4;
    
    %pressure field
    p = -(PX*fx+PY*fy)/pi/4 + visc2*(B1*ux+B2*uy)/pi/4;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    InOut = ones(numel(X0),1);
    u = PARAM.mass/2/pi/PARAM.alpha/PARAM.R;  %velocity in every point
    ux_inner = u.*r(1,n+m+m+j+1:end);
    uy_inner = u.*r(2,n+m+m+j+1:end);
    u_inner = zeros(2*PARAM.q,1);
    u_inner(1:2:end-1) = ux_inner;
    u_inner(2:2:end) = uy_inner;
    [~,~,U,V,p] = close_interface_visu...
        (Xsing(fixed+1:end),Ysing(fixed+1:end),X0,Y0,U,V,p,InOut,u_inner,2);

    %reshape in the grid shape
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    p = reshape(p,radial,axial);

    figure
    %subplot(2,1,1)
    quiver(X,Y,U,V)
    axis equal
    hold on
    %plot(a(1:fixed+1),b(1:fixed+1),'k',a(fixed+2:end),b(fixed+2:end),'k')
    plot(a(fixed+2:end),b(fixed+2:end),'k')
    %quiver([(a(1:fixed)+a(2:fixed+1))/2 (a(fixed+2:end-1)+a(fixed+3:end))/2],...
        %[(b(1:fixed)+b(2:fixed+1))/2 (b(fixed+2:end-1)+b(fixed+3:end))/2],ux',uy')
    if PARAM.continuity==1
        plot(PARAM.x_inj,PARAM.y_inj,'or')
    end
    hold off
    xlabel('x')
    ylabel('y')
    title(['flow field, \mu=' num2str(visc2)])
    
    figure
    %subplot(2,1,2)
    contourf(X,Y,p,200,'LineStyle','none')
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    title(['pressure field, \mu=' num2str(visc2)])
    hold on
    %plot(a(1:fixed+1),b(1:fixed+1),'k',a(fixed+2:end),b(fixed+2:end),'k')
    plot(a(fixed+2:end),b(fixed+2:end),'k')
    hold off
    
    
end