%micromotor evoultion till convergence (hopefully)

close all

%number of elements
q = PARAM.q;

%give reaches convergence
ind = find(risa(2,:)==0,1,'first');
DOplot = 0;

if isempty(ind)
    ind = PARAM.loop/PARAM.checkpoint;
end

%volume preallocation
V1 = zeros(ind,1);
V2 = zeros(ind,1);
xcm1 = zeros(ind-1,1);
xcm2 = zeros(ind-1,1);
dist = zeros(ind,1);
width1 = zeros(ind,1);
width2 = zeros(ind,1);

%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i = 1:ind-1
    
    disp([num2str(i) ' of ' num2str(ind-1)])
    
    if  sum(i==1:10:ind) && DOplot==1
    figure(1)
    plot(risa(1:q+1,i),risb(1:q+1,i),'or-',risa(q+2:end,i),risb(q+2:end,i),'ro-')
    hold on
    plot(risa(1:q+1,i),-risb(1:q+1,i),'or-',risa(q+2:end,i),-risb(q+2:end,i),'or-')
    %axis([min(a)-0.2 max(a)+0.2 -1 1]);
    title(['Ca1=' num2str(PARAM.Ca1) 'Ca2=' num2str(PARAM.Ca2) ' \lambda=' num2str(PARAM.visc)...
            ' \alpha1=' num2str(PARAM.alpha1) ' \alpha2=' num2str(PARAM.alpha2)])
    xlabel('x')
    ylabel('r')
    hold off
    axis equal
    legend('drop1','drop2')
    drawnow
    grid on
    end
    
    V1(i) = axis_int_gauss_vect(risa(1:PARAM.q+1,i)',risb(1:PARAM.q+1,i)');
    V2(i) = axis_int_gauss_vect(risa(PARAM.q+2:end,i)',risb(PARAM.q+2:end,i)');
    xcm1(i) = center_mass(risa(1:PARAM.q+1,i)',risb(1:PARAM.q+1,i)');
    xcm2(i) = center_mass(risa(PARAM.q+2:end,i)',risb(PARAM.q+2:end,i)');
    dist(i) = risa(PARAM.q+1,i)-risa(PARAM.q+2,i);
    width1(i) = max(risb(1:q+1,i));
    width2(i) = max(risb(q+2:end,i));
    
end

figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

%plot converged droplet shape
subplot(3,2,1)
plot(risa(1:PARAM.q+1,ind-1),risb(1:PARAM.q+1,ind-1),'b-',risa(PARAM.q+2:end,ind-1),risb(PARAM.q+2:end,ind-1),'r-')
hold on
axis([min(a)-0.2 max(a)+0.2 -1 1])
plot(risa(1:PARAM.q+1,ind-1),-risb(1:PARAM.q+1,ind-1),'b-',risa(PARAM.q+2:end,ind-1),-risb(PARAM.q+2:end,ind-1),'r-')
hold off
grid on
title(['Ca1=' num2str(PARAM.Ca1) 'Ca2=' num2str(PARAM.Ca2) ' \lambda=' num2str(PARAM.visc) ...
            ' \alpha1=' num2str(PARAM.alpha1) ' \alpha2=' num2str(PARAM.alpha2)])
xlabel('x')
ylabel('r')
hold off
legend('drop1','drop2')
axis equal

%plot velocity residuals
subplot(3,2,2)
%semilogy((2:ind-1)*PARAM.checkpoint,residuals(2:ind-1),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('res')
title('velocity and cm residuals')

%plot curvature in the final configuration
% subplot(2,2,3)
% plot(risa(1:PARAM.q+1,ind-1),KKK(1:PARAM.q+1,ind-2),'o-')
% hold on
% plot(risa(PARAM.q+2:end,ind-1),KKK(PARAM.q+2:end,ind-2),'o-')
% hold off
% title('curvature')
% xlabel('x')
% ylabel('k')
% grid on

%plot distnce between the two droplet
subplot(3,2,3)
plot((1:ind-3)*PARAM.checkpoint,dist(1:ind-3),'-','LineWidth',2)
xlabel('iteration')
ylabel('dist')
%title('Distance between interfaces on the axis')
grid on

%v_max poiseuille
v_max = 2*PARAM.Q/pi/PARAM.R^2;
v_cm1 = diff(xcm1(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;
v_cm2 = diff(xcm2(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;

%plot center of mass velocity
subplot(3,2,4)
%V_rel = abs(V-V(1))/V(1);
%figure
%semilogy((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
plot((1:ind-3)*PARAM.checkpoint,v_cm1,(1:ind-3)*PARAM.checkpoint,v_cm2,'-','LineWidth',2)
legend('drop1','drop2')
xlabel('iteration')
ylabel('vCM')
%title('Center of mass velocity')
grid on

%plot width of the two droplet
subplot(3,2,5)
plot((1:ind-3)*PARAM.checkpoint,width1(1:ind-3),(1:ind-3)*PARAM.checkpoint,width2(1:ind-3),'-','LineWidth',2)
xlabel('iteration')
legend('drop1','drop2')
ylabel('width')
%title('Width of the two droplets')
grid on

%plot volume error of the two droplet
V1ref = 4/3*pi*PARAM.R^3*PARAM.alpha1^3;
V2ref = 4/3*pi*PARAM.R^3*PARAM.alpha2^3;
subplot(3,2,6)
plot((1:ind-3)*PARAM.checkpoint,(V1(1:ind-3)-V1ref)/V1ref,(1:ind-3)*PARAM.checkpoint,(V2(1:ind-3)-V2ref)/V2ref,'-','LineWidth',2)
xlabel('iteration')
legend('drop1','drop2')
ylabel('errV')
%title('Width of the two droplets')
grid on

