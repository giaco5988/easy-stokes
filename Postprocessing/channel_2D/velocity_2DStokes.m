%velocity field for channel with mass injection trough potential flow

function [U,V,X,Y,p] = velocity_2DStokes(a,b,solution,resx,resy,PARAM)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

  n = PARAM.n;  m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
  
  wall = n+2*m+j;
  
  solution = solution(1:2*(n+2*m+j+q));
  
  computePress = 1;
  RelFrame = 1;

  visc2 = 1;
  visc1 = PARAM.visc;
  
  yy1 = max(b)-0.05;
  yy2 = min(b)+0.05;
  xx1 = min(a(n+2*m+j+2:end))-PARAM.R;
  xx2 = max(a(n+2*m+j+2:end))+PARAM.R;
  
%   yy1 = max(b)-0.5;
%   yy2 = min(b)+0.5;
%   xx1 = min(a)+4.5;
%   xx2 = max(a)-4.5;

    %interface normal and curvature
  [~, bx, cx, ~, ~, by, cy, ~] = my_spline_periodic ([a(n+m+j+m+2:end-1)' a(n+m+j+m+2)], [b(n+m+j+m+2:end-1)' b(n+m+j+m+2)]);
  
  %compute normals to the interface and curvature
  N = [by./sqrt(bx.*bx+by.*by); -bx./sqrt(bx.*bx+by.*by)];
  K = curv_spline_closed(bx,by,cx,cy);
  
  %stresses
  if PARAM.stress==1
    Q = 2/3*PARAM.press_grad/PARAM.visc*PARAM.R^3;
    Umax =  1.5*Q/2/PARAM.R;
  elseif PARAM.stress==0
    Umax = 1.5*PARAM.Q/2/PARAM.R;
  end
  gamma = PARAM.visc*Umax/PARAM.Ca;
  %gamma = PARAM.gamma;
  df = K*gamma;
  [fxDrop,fyDrop] = stress_diff(df,N,q);
  
  %mesh finess per unit length
  radial = round(resy*abs(yy1-yy2));
  axial = round(resx*abs(xx1-xx2));

    %create the grid
    x = linspace(xx1,xx2,axial);
    y = linspace(yy2,yy1,radial);

    [X,Y] = meshgrid(x,y);
    %X = 5;
    %Y = 0;
    
    X0 = X(:);
    Y0 = Y(:);
    
    %compute the necessary Green's function
    %[GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu_2DStokes(a,b,X0,Y0);
    [GXXconst,GXYconst,GYYconst,A11const,A12const,A22const] = GT_2DStokes_constElement(a(1:n+2*m+j+1),b(1:n+2*m+j+1),X0,Y0);
    [GXXlinSP,GXYlinSP,GYYlinSP,A11linSP,A12linSP,A22linSP] = GT_2DStokes_LinearSplinesElement(a(n+2*m+j+2:end),b(n+2*m+j+2:end),X0,Y0);
    
    %merge the two
    GXX = [GXXconst GXXlinSP];
    GXY = [GXYconst GXYlinSP];
    GYX = GXY;
    GYY = [GYYconst GYYlinSP];
    A11 = [A11const A11linSP];
    A12 = [A12const A12linSP];
    %A21 = A21;
    A22 = [A22const A22linSP];
       
    if computePress==1
        %in order to compute the pressure, INTEGRATIOMN ON WALL
        [PXwall,PYwall,PIXXwall,PIXYwall,PIYXwall,PIYYwall] = computeP_visu_2DStokes(a(1:wall+1),b(1:wall+1),X0,Y0);
        
        %in order to compute the pressure, INTEGRATIOMN ON WALL
        [PXdrop,PYdrop,PIXXdrop,PIXYdrop,PIYXdrop,PIYYdrop] = computeP_visu_2DStokesSpline(a(wall+2:end),b(wall+2:end),X0,Y0);
        
        PX = [PXwall PXdrop];
        PY = [PYwall PYdrop];
        PIXX = [PIXXwall PIXXdrop];
        PIXY = [PIXYwall PIXYdrop];
        PIYX = [PIYXwall PIYXdrop];
        PIYY = [PIYYwall PIYYdrop];
    end

    %add velocity field and stresses with potential flow
    u = 0;
    v = 0;
    ux = 0;
    uy = 0;
    NXX = 0;
    NXY = 0;
    NYX = 0;
    NYY = 0;
    for l = 1:numel(PARAM.mass)

            %in order to modify the physical boundary conditions
            [~,u_temp,v_temp,nxx_temp,nxy_temp,nyx_temp,nyy_temp] = gf_2D_Laplace((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,...
                PARAM.x_inj(l),PARAM.y_inj(l),PARAM.R,-PARAM.R,visc2,0,PARAM.mass(l));

            ux = ux + u_temp;
            uy = uy + v_temp;
            NXX = NXX + nxx_temp;
            NXY = NXY + nxy_temp;
            NYX = NYX + nyx_temp;
            NYY = NYY + nyy_temp;

            %compute the velocities from potential flow that I will overlap
            [~,u_temp,v_temp] = gf_2D_Laplace(X0,Y0,...
                PARAM.x_inj(l),PARAM.y_inj(l),PARAM.R,-PARAM.R,visc2,PARAM.p_out,PARAM.mass(l));

            u = u + u_temp;
            v = v + v_temp;

    end
    
    %normal vector
    no = [sqrt(diff(b(1:wall+1)).^2+diff(a(1:wall+1)).^2); sqrt(diff(b(wall+2:end)).^2+diff(a(wall+2:end)).^2)];
    r = [-diff(b(1:2*m+n+j+1))'./no(1:wall)' -diff(b(2*m+n+j+2:end))'./no(wall+1:end)';  diff(a(1:2*m+n+j+1))'./no(1:wall)' diff(a(2*m+n+j+2:end))'./no(wall+1:end)']';
    
    if computePress==1
        R1 = repmat(r(1:end,1)',numel(X0),1);
        R2 = repmat(r(1:end,2)',numel(X0),1);
    end
    
    %double layer
%     A11 = TXXX.*R1+TXXY.*R2;
%     A12 = TXYX.*R1+TXYY.*R2;
%     A21 = TYXX.*R1+TYXY.*R2;
%     A22 = TYYX.*R1+TYYY.*R2;
    
    if computePress==1
        %double layer pressure
        B1 = PIXX.*R1+PIXY.*R2;
        B2 = PIYX.*R1+PIYY.*R2;
    end
    
    %stresses due to potential flow
    rPot = [-diff(b)./sqrt(diff(a).^2+diff(b).^2) diff(a)./sqrt(diff(a).^2+diff(b).^2)];
    fx = NXX.*rPot(:,1)+NXY.*rPot(:,2);
    fy = NYX.*rPot(:,1)+NYY.*rPot(:,2);
    
    %BC inlet
    Ysing = -diff(b(n+m+1:n+m+j+1))/2;
    uIN = 3/4*PARAM.Q/PARAM.R*(1-Ysing.^2/PARAM.R^2);
    
    %MODIFIED BOUNDARY CONDITIONS AND COMPLEMENTARY SOLUTION (from BEM)
    if PARAM.continuity==1
        uux = [solution(1:2:2*n-1)-ux(1:n)'; -ux(n+1:m+n)'; solution(2*(n+m)+1:2:2*(n+m+j)-1)-ux(n+m+1:n+m+j)'; -ux(n+m+j+1:end)'];
        uuy = zeros(2*m+n+j,1);
    else
        if PARAM.stress==0
        uux = [solution(1:2:2*n-1); zeros(m,1); uIN; zeros(m,1); solution(2*(n+2*m+j)+1:2:end)];
        uuy = [zeros(2*m+n+j,1); solution(2*(n+2*m+j)+2:2:end)];
        end
    end
    
    if PARAM.continuity==1
        ffx = [PARAM.p_out-fx(1:n)'; solution(2*n+1:2:2*(n+m)-1,1)-fx(n+1:n+m)'; ...
            -(PARAM.press_grad+PARAM.p_out)-fx(n+m+1:n+m+j)'; solution(2*(n+m+j)+1:2:end-1)-fx(n+m+j+1:end)'];
        ffy = solution(2:2:end)-fy';
    else
        if PARAM.stress==0
        ffx = [PARAM.p_out*ones(n,1); solution(2*n+1:2:2*(n+2*m+j)-1,1); fxDrop'];
        ffy = [solution(2:2:2*(n+2*m+j)); fyDrop'];
        end
    end
    
    %flow field from complementary solution
    U = (-GXX*ffx-GXY*ffy + visc2*(A11(:,1:n+2*m+j)*uux(1:n+2*m+j)+A12(:,1:n+2*m+j)*uuy(1:n+2*m+j)) + (visc2-visc1)*(A11(:,n+2*m+j+1:end)*uux(n+2*m+j+1:end)+A12(:,n+2*m+j+1:end)*uuy(n+2*m+j+1:end)))/pi/4;
    V = (-GYX*ffx-GYY*ffy + visc2*(A12(:,1:n+2*m+j)*uux(1:n+2*m+j)+A22(:,1:n+2*m+j)*uuy(1:n+2*m+j)) + (visc2-visc1)*(A12(:,n+2*m+j+1:end)*uux(n+2*m+j+1:end)+A22(:,n+2*m+j+1:end)*uuy(n+2*m+j+1:end)))/pi/4;
    
    %pressure field
    if computePress==1
        %FAKE 0th order velocity
        %UdropConst = (uux(n+2*m+j+1:end)+[uux(n+2*m+j+2:end); uux(1)])/2;
        %VdropConst = (uuy(n+2*m+j+1:end)+[uuy(n+2*m+j+2:end); uuy(1)])/2;
        p = -(PX*ffx+PY*ffy)/pi/4 + (visc2*(B1(:,1:n+2*m+j)*uux(1:n+2*m+j)+B2(:,1:n+2*m+j)*uuy(1:n+2*m+j)) + (visc2-visc1)*(B1(:,n+2*m+j+1:end)*uux(n+2*m+j+1:end)+B2(:,n+2*m+j+1:end)*uuy(n+2*m+j+1:end)))/pi/4;
    end
    
    %physical flow field adding the potential
    if PARAM.continuity==1
        U = U + u;
        V = V + v;
    end
    
    %check if evaluation point is very close to the interface and change it if
    %it is so switch to solution data
    %p = 1;
    InOut = ones(numel(X0),1);
    [X0,Y0,~,~,p,InOut] = close_interface_visu(a(wall+2:end),b(wall+2:end),X0,Y0,U,V,p,InOut,solution(2*wall+1:end),4);

    %reshape in the grid shape
    X = reshape(X0,radial,axial);
    Y = reshape(Y0,radial,axial);
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    if computePress==1
        p = reshape(p,radial,axial);
    end

    figure
    xDrop = a(n+2*m+j+2:end);   yDrop = b(n+2*m+j+2:end);
    xWall = a(2:n+2*m+j+1); yWall = b(2:n+2*m+j+1);
    %subplot(2,1,1)
    if RelFrame==1
        U = U-uux(n+2*m+j+1);
    end
    quiver(X,Y,U,V)
    axis equal
    axis([min(a(n+2*m+j+2:end))-PARAM.R max(a(n+2*m+j+2:end))+PARAM.R min(b)-0.1 max(b)+0.1])
    hold on
    plot(xWall,yWall,'k','LineWidth',2)
    plot(xDrop,yDrop,'k','LineWidth',2)
    %quiver((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,ux+uux',uy+uuy')
    if PARAM.continuity==1
        plot(PARAM.x_inj,PARAM.y_inj,'or')
    end
    %plot(a,b,'k')
    hold off
    xlabel('x')
    ylabel('y')
    title(['flow field, \lambda=' num2str(visc2) ' Ca=' num2str(PARAM.Ca)])
    
%     figure
%     %plot velocity magnitude
%     subplot(2,1,1)
%     contourf(X,Y,sqrt(U.*U+V.*V),200,'LineStyle','none')
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     title(['velocity magnitude, \mu=' num2str(visc2)])
%     hold on
%     plot(a,b)
%     hold off
    
    if computePress==1
    %plot pressure
    %subplot(2,1,2)
    figure
    contourf(X,Y,p,200,'LineStyle','none')
    colorbar
    hold on
    plot(xWall,yWall,'k','LineWidth',2)
    plot(xDrop,yDrop,'k','LineWidth',2)
    axis equal
    axis([min(a(n+2*m+j+2:end))-PARAM.R max(a(n+2*m+j+2:end))+PARAM.R min(b)-0.1 max(b)+0.1])
    xlabel('x')
    ylabel('y')
    title(['pressure field, \lambda=' num2str(visc2) ' Ca=' num2str(PARAM.Ca)])
    hold off
    end
    
end