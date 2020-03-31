%micromotor evoultion till convergence (hopefully)

close all

set(0,'defaultaxesfontsize',15,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

%number of elements
n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;    p = PARAM.p;

%options
plotRT = 0;

%give reaches convergence
ind = find(risa(2,:)==0,1,'first');

if isempty(ind)==1
    ind = PARAM.loop/PARAM.checkpoint;
end

%volume preallocation
V1 = zeros(ind,1);
V2 = zeros(ind,1);
v1Front = zeros(ind,1);
v2Front = zeros(ind,1);
xcm1 = zeros(ind,1);
xcm2 = zeros(ind,1);
width1 = zeros(ind,1);
width2 = zeros(ind,1);

%v_max poiseuille
v_max = 2*PARAM.Q/pi/PARAM.R^2;

r1 = PARAM.R*PARAM.alpha1;
r2 = PARAM.R*PARAM.alpha2;

%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i = 1:ind-1
    disp(i)
    
    %wall nodes
    aWall = risa(1:n+m+j+1,i);  bWall = risb(1:n+m+j+1,i);
    
    %drops nodes
    aDrop1 = risa(n+m+j+2:n+m+j+q+2,i); bDrop1 = risb(n+m+j+2:n+m+j+q+2,i);
    aDrop2 = risa(n+m+j+q+3:n+m+j+q+p+3,i); bDrop2 = risb(n+m+j+q+3:n+m+j+q+p+3,i);
    
    if plotRT==1
    figure(1)
    plot(aWall,bWall,'ob-',aWall,-bWall,'ob-')
    hold on
    plot(aDrop1,bDrop1,'ro-',aDrop1,-bDrop1,'ro-')
    plot(aDrop2,bDrop2,'ro-',aDrop2,-bDrop2,'ro-')
    title(['Ca1=' num2str(PARAM.Ca1) ' \lambda_1=' num2str(PARAM.visc) ' \alpha_1=' num2str(PARAM.alpha1)...
        'Ca2=' num2str(PARAM.Ca2) ' \lambda_2=' num2str(PARAM.visc) ' \alpha_2=' num2str(PARAM.alpha2)])
    xlabel('x')
    ylabel('r')
    hold off
    grid on
    if PARAM.dist>0
    axis([PARAM.L/2-5*r2 PARAM.L/2+5*r1 -1 1])
    else
    axis([PARAM.L/2-2*r1 PARAM.L/2+2*r1 -1 1])
    end
    axis equal
    drawnow
    hold off
    end
    
    V1(i) = axis_int_gauss_vect(aDrop1',bDrop1');
    V2(i) = axis_int_gauss_vect(aDrop2',bDrop2');
    v1Front(i) = risy(2*(n+m+j)+1,i)/v_max;
    v2Front(i) = risy(2*(n+m+j+q+1)+1,i)/v_max;
    xcm1(i) = center_mass(aDrop1,bDrop1);
    xcm2(i) = center_mass(aDrop2,bDrop2);
    width1(i) = 2*max(bDrop1);
    width2(i) = 2*max(bDrop2);
    
end

figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

%plot converged droplet shape
subplot(3,2,1)
plot(aWall,bWall,'ob-',aWall,-bWall,'ob-')
hold on
plot(aDrop1,bDrop1,'ro-',aDrop1,-bDrop1,'ro-')
plot(aDrop2,bDrop2,'ro-',aDrop2,-bDrop2,'ro-')
title(['Ca1=' num2str(PARAM.Ca1) ' \lambda_1=' num2str(PARAM.visc) ' \alpha_1=' num2str(PARAM.alpha1)...
        'Ca2=' num2str(PARAM.Ca2) ' \lambda_2=' num2str(PARAM.visc) ' \alpha_2=' num2str(PARAM.alpha2)])
xlabel('x')
ylabel('r')
hold off
grid on
if PARAM.dist>0
    axis([PARAM.L/2-5*r2 PARAM.L/2+5*r1 -1 1])
else
    axis([PARAM.L/2-2*r1 PARAM.L/2+2*r1 -1 1])
end
axis equal


%v_cm1 = diff(xcm1(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;
%v_cm2 = diff(xcm2(1:end-1))/PARAM.deltaT/PARAM.checkpoint/v_max;

%plot center of mass velocity
subplot(3,2,4)
plot((1:ind-2)*PARAM.checkpoint,v1Front(1:ind-2),(1:ind-2)*PARAM.checkpoint,v2Front(1:ind-2),'-','LineWidth',2)
legend('drop1','drop2')
xlabel('iteration')
ylabel('vCM')
grid on

%plot velocity residuals
subplot(3,2,2)
try
    semilogy((2:ind-1)*PARAM.checkpoint,residuals(2:ind-1),'LineWidth',2)
catch exception
    semilogy((2:ind-1)*PARAM.checkpoint,many_res(2:ind-1),'LineWidth',2)
end
grid on
xlabel('iteration')
ylabel('res')
title('velocity residuals')

%plot distance between center of mass
subplot(3,2,3)
loglog((1:ind-1)*PARAM.checkpoint,xcm1(1:ind-1)-xcm2(1:ind-1),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('dist x_{cm}')
%title('velocity residuals')

%plot curvature in the final configuration
% subplot(2,2,3)
% plot(risa(n+m+j+2:n+m+j+q+2,ind-1),KKK(1:q+1,ind-2),'o-b',risa(n+m+j+q+3:n+m+j+q+p+3,ind-1),KKK(q+2:q+p+2,ind-2),'o-b')
% title('curvature')
% xlabel('x')
% ylabel('k')
% grid on

%plot width of the two droplet
subplot(3,2,5)
plot((1:ind-3)*PARAM.checkpoint,width1(1:ind-3),(1:ind-3)*PARAM.checkpoint,width2(1:ind-3),'-','LineWidth',2)
xlabel('iteration')
legend('drop1','drop2')
ylabel('width')
%title('Width of the two droplets')
grid on

subplot(3,2,6)
%plot volume variation
V_rel1 = abs(V1-V1(1))/V1(1);
V_rel2 = abs(V2-V2(1))/V2(1);
%figure
%semilogy((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
plot((2:ind-1)*PARAM.checkpoint,V_rel1(2:ind-1),'LineWidth',2)
hold on
plot((2:ind-1)*PARAM.checkpoint,V_rel1(2:ind-1),'LineWidth',2)
hold off
legend('drop1','drop2')
xlabel('iteration')
ylabel('errV_{rel}')
title('Volume variation')
grid on


