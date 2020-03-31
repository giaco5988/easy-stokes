%find stagnation point along the axis

%close all

%number of loops
finish = find(risb(2,:)==0,1,'first')-2;    %stop before crash
if isempty(finish)
    finish = loop/checkpoint;
end

%stop at this loop
stop = 100;

%allocation
dist = zeros(finish-2,1);

for i = 3:finish
    
    %i = 120;
    
    disp(i-2)
    
    %i = 10;

    %solution at this step
    solution = risy(:,i);
    
    %interface shape at this iteration
    a_bef = risa(:,i-2)';
    b_bef = risb(:,i-2)';
    a = risa(:,i-1)';
    b = risb(:,i-1)';
    
    %smallest element size
    all_ds = sqrt((a(1:end-1)-a(2:end)).^2+(b(1:end-1)-b(2:end)).^2);
    %ds = min(all_ds);
    %first element
    ds = all_ds(1);
    
    %center of mass and its velocity
    xcm_old = center_mass(a_bef,b_bef);
    xcm = center_mass(a,b);
    v_centermass = -(xcm-xcm_old)/deltaT/checkpoint;
    a = a-xcm;
    
    %compute the spline coeff
    [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
    
    %BC stresses
    [df_x,df_y] = df_buoyancy_spline_symm(a,b,Ca,visc,1);
    
    %BC velocities
    ux_sol = solution(1:2:end-1);
    uy_sol = solution(2:2:end);
    
    space = ds/2;     %spacing close to the interface
    %n = 5; %number of points
    %t = [0:0.1:2]';
    %t = [linspace(0,risa(1,i),n) linspace(risa(1,i),2,n)]'; t = [t(1:n); t(n+2:end)];
    t = [0:0.2:0.9*a(1) 0.9*a(1):space:a(1)-space a(1)+space:space:1.1*a(1) 1.1*a(1):0.2:2*a(1)]';
    %t1 = [0:0.2:a(1)-20*ds]';  t2 = [a(1)-20*ds:ds:a(1)+20*ds]'; t3 = [a(1)+20*ds:ds:2*a(1)]';
    X0 = t; %constant spacing
    %X1 = t1;    X2 = t2;    X3 = t3; %denser spacing around the interface
    %X0 = [X1; X2; X3];
    Y0 = zeros(numel(X0),1);
    %Y1 = zeros(numel(X1),1);
    %Y2 = zeros(numel(X2),1);
    %Y3 = zeros(numel(X3),1);
    
    %figure out if I'm inside or outside the droplet
    InOut = FigInOut(a,b,X0,Y0);
    InOut = (InOut>pi)/visc + (InOut<pi);

    %compute the necessary Green's function
    [GXX,GXY,GYX,GYY,A11,A12,A21,A22] = computeGT_spline_visu(ax,ay,bx,by,cx,cy,dx,dy,X0',Y0');
    
    %compute the integral to obtain the velocities
    U = (-GXX*df_x' - GXY*df_y' + (1-visc)*A11*ux_sol + (1-visc)*A12*uy_sol)/8/pi;
    %V = (-GYX*df_x' - GYY*df_y' + (1-visc)*A21*ux_sol + (1-visc)*A22*uy_sol)/8/pi;
    
    %check if evaluation point is very close to the interface and change it if
    %it is so
    p = 0;
    [X0,Y0,U,V,p,InOut] = close_interface_visu(a,b,X0,Y0,U,V,p,InOut,solution,10);
    
    %moltiply time 1/lambda the points inside the droplet
    U = InOut.*U;
    %V = InOut.*V;
    
    %stagnation point
    %x = linspace(X0(1),X0(end),10000); %finer grid for interpolation
    [~,ind1] = min(abs(X0-a(1)));
    x1 = linspace(X0(1),X0(ind1),10000); %finer grid for interpolation
    x2 = linspace(X0(ind1+1),X0(end),10000); %finer grid for interpolation
    %u_spline = spline(X0,U,x);
    u_spline1 = spline(X0(1:ind1),U(1:ind1),x1);
    u_spline2 = spline(X0(ind1+1:end),U(ind1+1:end),x2);
    x = [x1 x2];
    u_spline = [u_spline1 u_spline2];
    [~,ind2] = min(abs(-u_spline-v_centermass));
    stagnation = x(ind2);
    
    %distance between rear tip and stagnation point
    dist(i-2) = a(1)-stagnation;
    
    figure(1)
    subplot(2,1,1)
    plot(b,-a,'b',-b,-a,'b','LineWidth',2)
    hold on
    axis([-0.2 0.2 -stagnation-0.2 -stagnation+0.2])
    %plot(Y0,-X0,'or')
    plot(0,-stagnation,'or')
    hold off
    axis equal
    xlabel('r')
    ylabel('x')
    title(strcat('t=',num2str(i*deltaT*checkpoint)))
    
    %figure(2)
    subplot(2,1,2)
    plot(-x,-u_spline-v_centermass,'LineWidth',2)
    %find tip
    [~,ind3] = min(abs(a(1)-x));
    hold on
    plot(-x(ind3),-u_spline(ind3)-v_centermass,'or')
    hold off
    xlabel('x')
    ylabel('vel')
    grid on

    if i==stop
        break;
    end
    
    %plot rear tip velocity (check if it oscillates)
%     figure(3)
%     hold on
%     plot(i,risy(1,i),'o')
%     grid on
%     hold off
  
end

%physical time
time = deltaT*checkpoint*2:deltaT*checkpoint:deltaT*checkpoint*(finish-1);

if stop>i
    stop = i;
end

figure
plot(time(1:stop-3),dist(1:stop-3),'LineWidth',2)
%plot(time,dist,'ro')
%semilogy(dist)
xlabel('time')
ylabel('dist')
grid on
title(['Stagnation point \lambda=' num2str(visc) ', Ca=' num2str(Ca) ', \delta=' num2str(D)])
%title(['Stagnation point \lambda=' num2str(visc) ', Ca=' num2str(Ca)])