%micromotor evoultion till convergence (hopefully)

close all

set(0,'defaultaxesfontsize',14,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%number of elements
q = PARAM.q;

%give reaches convergence
ind = find(risa(2,:)==0,1,'first');

if isempty(ind)
    ind = PARAM.loop/PARAM.checkpoint;
end

%volume preallocation
V = zeros(ind-1,1);
xcm = zeros(ind-1,1);   %-2 because I have the break of the convergenge
width = zeros(ind-1,1);

%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i = 1:ind-1
    disp(i)
    
    if sum(1:1:ind==i)
    figure(1)
    plot(risa(1:q+1,i),risb(1:q+1,i),'ro-',risa(1:q+1,i),-risb(1:q+1,i),'ro-')
    title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc)...
            ' \alpha=' num2str(PARAM.alpha)])
    xlabel('x')
    ylabel('r')
    hold off
    axis equal
    end
    
    V(i) = axis_int_gauss_vect(risa(:,i)',risb(:,i)');
    xcm(i) = center_mass(risa(:,i)',risb(:,i)');
    width(i) = max(risb(1:q+1,i));
    
end

%V = axis_int_gauss_vect(risa',risb');
%xcm = center_mass(risa',risb');

figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

%plot converged droplet shape
subplot(3,2,1)
plot(risa(1:PARAM.q+1,ind-1),risb(1:PARAM.q+1,ind-1),'ro-',risa(1:PARAM.q+1,ind-1),-risb(1:PARAM.q+1,ind-1),'ro-')
hold on
plot(risa(PARAM.q+2:end,ind-1),risb(PARAM.q+2:end,ind-1),'ro-',risa(PARAM.q+2:end,ind-1),-risb(PARAM.q+2:end,ind-1),'ro-')
hold off
title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ...
            ' \alpha=' num2str(PARAM.alpha)])
xlabel('x')
ylabel('r')
hold off
axis equal

%plot velocity residuals
subplot(3,2,2)
semilogy((2:ind-1)*PARAM.checkpoint,residuals(2:ind-1),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('res')
title('velocity residuals')

%plot width in the final configuration
subplot(3,2,3)
%plot(risa(1:PARAM.q+1,ind-1),KKK(1:PARAM.q+1,ind-2),'o-')
plot((1:ind-1)*PARAM.checkpoint,width)
hold on
%plot(risa(PARAM.q+2:end,ind-1),KKK(PARAM.q+2:end,ind-2),'o-')
hold off
title('droplet width')
xlabel('x')
ylabel('k')
grid on

%v_max poiseuille
v_max = 2*PARAM.Q/pi/PARAM.R^2;
v_cm = diff(xcm(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;

subplot(3,2,4)
plot((1:ind-3)*PARAM.checkpoint,v_cm,'-','LineWidth',2)
xlabel('iteration')
ylabel('v')
title('Center of mass velocity/v_max')
grid on

%plot volume error of the two droplet
subplot(3,2,5)
plot((1:ind-3)*PARAM.checkpoint,V(1:ind-3),'-','LineWidth',2)
xlabel('iteration')
ylabel('V')
title('volume')
grid on

% subplot(2,2,4)
% %plot volume variation
% V_rel = abs(V-V(1))/V(1);
% %figure
% %semilogy((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
% plot((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
% xlabel('iteration')
% ylabel('V_{rel}')
% title('Volume variation')
% grid on