%velocity field for channel with mass injection trough potential flow

function [U,V,p,X,Y] = velocity_2DStokes(a,b,solution,resx,resy,PARAM)

  set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

  n = PARAM.n;  m = PARAM.m;    j = PARAM.j;

  visc2 = PARAM.visc2;
  
  yy1 = max(b)-0.1;
  yy2 = min(b)+0.1;
  xx1 = min(a)+0.1;
  xx2 = max(a)-0.1;
  
%   yy1 = max(b)-0.5;
%   yy2 = min(b)+0.5;
%   xx1 = min(a)+4.5;
%   xx2 = max(a)-4.5;
  
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
    [GXX,GXY,GYX,GYY,TXXX,TXXY,TXYX,TXYY,TYXX,TYXY,TYYX,TYYY] = computeGT_visu_2DStokes(a,b,X0,Y0);
       
    %in order to compute the pressure
    [PX,PY,PIXX,PIXY,PIYX,PIYY] = computeP_visu_2DStokes(a,b,X0,Y0);

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
    no = sqrt(diff(b).^2+diff(a).^2);
    r = [-diff(b)./no; diff(a)./no];
    
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
    
    %stresses due to potential flow
    fx = NXX.*r(1,:)+NXY.*r(2,:);
    fy = NYX.*r(1,:)+NYY.*r(2,:);
    
    %MODIFIED BOUNDARY CONDITIONS AND COMPLEMENTARY SOLUTION (from BEM)
    uux = [solution(1:2:2*n-1)-ux(1:n)'; -ux(n+1:m+n)'; solution(2*(n+m)+1:2:2*(n+m+j)-1)-ux(n+m+1:n+m+j)'; -ux(n+m+j+1:end)'];
    uuy = zeros(2*m+n+j,1);
    
    ffx = [PARAM.p_out-fx(1:n)'; solution(2*n+1:2:2*(n+m)-1,1)-fx(n+1:n+m)'; ...
        -(PARAM.press_grad+PARAM.p_out)-fx(n+m+1:n+m+j)'; solution(2*(n+m+j)+1:2:end-1)-fx(n+m+j+1:end)'];
    ffy = solution(2:2:end)-fy';
    
    %flow field from complmentary solution
    U = (-GXX*ffx-GXY*ffy + visc2*(A11*uux+A12*uuy))/pi/4;
    V = (-GYX*ffx-GYY*ffy + visc2*(A21*uux+A22*uuy))/pi/4;
    
    %pressure field
    p = -(PX*ffx+PY*ffy)/pi/4 + (visc2*(B1*uux+B2*uuy))/pi/4;
    
    %physical flow field adding the potential
    if PARAM.continuity==1
        U = U + u;
        V = V + v;
    end

    %reshape in the grid shape
    U = reshape(U,radial,axial);
    V = reshape(V,radial,axial);
    p = reshape(p,radial,axial);

    figure
    %subplot(2,1,1)
    quiver(X,Y,U,V)
    axis equal
    hold on
    %quiver((a(1:end-1)+a(2:end))/2,(b(1:end-1)+b(2:end))/2,ux+uux',uy+uuy')
    if PARAM.continuity==1
        plot(PARAM.x_inj,PARAM.y_inj,'or')
    end
    %plot(a,b,'k')
    hold off
    xlabel('x')
    ylabel('y')
    title(['flow field, \mu=' num2str(visc2)])
    
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
    
    %plot pressure
    %subplot(2,1,2)
    figure
    contourf(X,Y,p,200,'LineStyle','none')
    colorbar
    axis equal
    xlabel('x')
    ylabel('y')
    title(['pressure field, \mu=' num2str(visc2)])
    hold on
    %plot(a,b,'k')
    hold off
    
end