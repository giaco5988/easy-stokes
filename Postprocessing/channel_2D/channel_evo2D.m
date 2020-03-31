%micromotor evoultion till convergence (hopefully)

close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%number of elements
n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;

%options
plotRT = 0;
stop = 0;

%give reaches convergence
ind = find(risa(2,:)==0,1,'first');
if isempty(ind)
    ind = PARAM.loop/PARAM.checkpoint;
end

if stop==1
br = 450;
ind = br;
end

%volume preallocation
V = zeros(ind,1);
velDrop = zeros(ind,1);
constH = zeros(ind,1);
minH = zeros(ind,1);

%Uavg of the inlet
Uavg = PARAM.Q/2/PARAM.R;

%figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i = 1:ind-1
    disp(i)
    
    %find wall
    indWall = find(abs((risb(:,i)-risb(1,i)).^2+(risa(:,i)-risa(1,i)).^2)<1e-6,2,'first');
    aWall = risa(1:indWall(2),i);  bWall = risb(1:indWall(2),i);
    
    %find drop
    indDrop = find(abs((risb(:,i)-risb(indWall(2)+1,i)).^2+(risa(:,i)-risa(indWall(2)+1,i)).^2)<1e-6,2,'first');
    aDrop = risa(indWall(2)+1:indDrop(2),i);   bDrop = risb(indWall(2)+1:indDrop(2),i);
       
    if plotRT==1
    figure(1)
    plot(aWall,bWall,'bo-')
    hold on
    plot(aDrop,bDrop,'ro-')
    title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ' \alpha=' num2str(PARAM.alpha)])
    xlabel('x')
    ylabel('r')
    hold off
    grid on
    hold off
    axis([-3*PARAM.R*PARAM.alpha+PARAM.L/2 3*PARAM.R*PARAM.alpha+PARAM.L/2 -1 1])
    axis equal
    drawnow
    end
    
    if max(bDrop)>PARAM.R
        warning('Interface escaped the domain')
        break;
    end
    
    V(i) = compute_area_2D(aDrop,bDrop);
    velDrop(i) = risy(2*n+4*m+2*j+1,i)/Uavg;
    constH(i) = min(max(risb(:,ind-1))-bDrop(q/4))/2/PARAM.R;
    minH(i) = min(max(risb(:,ind-1))-bDrop)/2/PARAM.R;
    
end

figure
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

%plot converged droplet shape
subplot(2,2,1)
plot(aWall,bWall,'ob-')
hold on
plot(aDrop,bDrop,'ro-')
title(['Ca=' num2str(PARAM.Ca) ' \lambda=' num2str(PARAM.visc) ' \alpha=' num2str(PARAM.alpha)])
xlabel('x')
ylabel('r')
grid on
hold off
axis([-3*PARAM.R*PARAM.alpha+PARAM.L/2 3*PARAM.R*PARAM.alpha+PARAM.L/2 -1 1])
axis equal

%plot velocity residuals
% subplot(2,2,2)
% try
%     semilogy((2:ind-1)*PARAM.checkpoint,residuals(2:ind-1),'LineWidth',2)
% catch exception
%     %semilogy((2:ind-1)*PARAM.checkpoint,many_res(2:ind-1),'LineWidth',2)
%     warning('No data for residuals')
% end
% grid on
% xlabel('iteration')
% ylabel('res')
% title('velocity residuals')

%plot droplet velocity
subplot(2,2,2)
semilogx((2:ind-2)*PARAM.checkpoint,velDrop(2:ind-2),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('v_{drop}')
title('drop velocity')

%plot film thickness
subplot(2,2,3)
semilogx((2:ind-2)*PARAM.checkpoint,constH(2:ind-2),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('h')
title('film thickness')

%plot film thickness
subplot(2,2,4)
semilogx((2:ind-2)*PARAM.checkpoint,minH(2:ind-2),'LineWidth',2)
grid on
xlabel('iteration')
ylabel('h')
title(' minimum film thickness')

%plot curvature in the final configuration
% subplot(2,3,5)
% coord = sqrt(diff(risa(n+2*m+j+2:end,ind-1)).^2 + diff(risb(n+2*m+j+2:end,ind-1)).^2);
% coord = [0; cumsum(coord)];
% plot(coord,[KKK(1:q,ind-2); KKK(1,ind-2)],'o-')
% title('curvature')
% xlabel('x')
% ylabel('k')
% grid on
% 
% subplot(2,3,6)
% %plot volume variation
% V_rel = abs(V-V(1))/V(1);
% %figure
% %semilogy((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
% plot((2:ind-1)*PARAM.checkpoint,V_rel(2:ind-1),'LineWidth',2)
% xlabel('iteration')
% ylabel('V_{rel}')
% title('Volume variation')
% grid on