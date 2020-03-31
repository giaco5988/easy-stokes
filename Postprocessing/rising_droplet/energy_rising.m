%compare kinetic and surface energy

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%load data
load('~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/prolate/q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.068_Ca=6_RK2.mat')

stop = 300;
step = 20;
plotNOW = 0;
plotVEL = 0;
plotDIFF = 0;
plotShapes = 0;

%geometry at first
[aIN,bIN] = draw_circle_lean(0,0,q,D);

%time at which plot the droplet shape
ttt = [3 6 10];

%plot shapes for those time
if plotShapes==1
figure
plot([aIN flip(aIN)],[bIN -flip(bIN)],'k')
hold on
for k = 1:numel(ttt)
   
    ite = ttt(k)/deltaT/checkpoint;
    a = risa(:,ite)';
    b = risb(:,ite)';
    
    plot([a flip(a)],[b -flip(b)])
end
axis equal
hold off
end

%for vel field
res = 40;    xlim = 2.5;  ylim = 2.5;
x = linspace(-xlim,xlim,2*res); y = linspace(0,ylim,res);
[X,Y] = meshgrid(x,y);
dx = diff(x);   dy = diff(Y);

%plot first domain
% figure(1)
% plot(risa(:,1),risb(:,1),'b',risa(:,1),-risb(:,1),'b')
% hold on
% plot(X,Y,'.r')
% plot(X,-Y,'.r')
% axis equal
% grid on
% axis([-xlim xlim -ylim ylim])
% xlabel('x')
% ylabel('r')
% hold off
% drawnow

Vin = axis_int_gauss(risa(:,1)',risb(:,1)');
Rin = nthroot(Vin/4*3/pi,3);

%initialize
Area = zeros(stop+1,1);
Area(1) = surf_gauss_vect(aIN,bIN);
magn = zeros(stop,1);
visc_diss = zeros(stop,1);
for i = 1:stop
    
    disp([num2str(i) ' of ' num2str(stop)])
    
    a = risa(:,i)';
    b = risb(:,i)';
    
    %plot droplet
    if plotNOW==1 && plotVEL==0
        figure(1)
        plot(a,b,'b',a,-b,'b')
        axis equal
        grid on
        axis([-xlim xlim -ylim ylim])
        xlabel('x')
        ylabel('r')
        drawnow
    end
    
    %compute area
    Area(i+1) = surf_gauss_vect(a,b);
    
    %compute vel field and magnitude
%     [U,V] = velocity_field_relaxation_tracers(a,b,risy(:,i),visc,Ca,X,Y);
%     both = U.^2+V.^2;
%     bothUP = both(1:end-1,:);   bothDOWN = both(2:end,:);   bothLEFT = both(:,1:end-1); bothRIGHT = both(:,2:end);
%     IntY = sum((bothUP.*Y(1:end-1,:)+bothDOWN.*Y(2:end,:)).*dy)/2;
%     magn(i) = sum((IntY(1:end-1)+IntY(2:end)).*dx)/2;
    
    %plot drolet and velocity field
    if plotVEL==1
        figure(1)
        plot(a,b,'b',a,-b,'b')
        hold on
        quiver(X,Y,U,V,'b')
        quiver(X,-Y,U,-V,'b')
        axis equal
        grid on
        axis([-xlim xlim -ylim ylim])
        xlabel('x')
        ylabel('r')
        hold off
        drawnow
    end
    
    %compute viscous dissipation on the surface
    ux = risy(1:2:end-1,i)';  uy = risy(2:2:end,i)';    %velocities
    [df_x,df_y] = df_surf_tens_spline_symm(a,b,Ca);  %BC stresses
    int = df_x.*ux + df_y.*uy;
    int1 = int(1:end-1);    int2 = int(2:end);  dl = sqrt(diff(b).^2+diff(a).^2);
    INT = pi*(int1.*b(1:end-1)+int2.*b(2:end)).*dl;
    visc_diss(i) = sum(INT);
    
%     figure(2)
%     plot(ux)
%     hold on
%     plot(uy)
%     hold off
%     xlabel('x')
%     ylabel('r')
%     grid on
%     drawnow
    
end

%compute energy
gamma = 1;
Es = Area-4*pi*Rin^2;   %surface energy
Ek = magn;   %kinetic energy
Pdiss = visc_diss;
Ps = -diff(Es)/deltaT;

time = deltaT*checkpoint:deltaT*checkpoint:deltaT*checkpoint*stop;

figure
plot(time,Ps)
hold on
plot([time(1:step:end) time(end)],-[Pdiss(1:step:end); Pdiss(end)],'o')
hold off
grid on
xlabel('t')
ylabel('Power')
legend('dS/st','Dissipation')

%plt dots for corresponding shapes
hold on
plot(0,-Pdiss(1),'.k','MarkerSize',50)
for k = 1:numel(ttt)
   
    ite = ttt(k)/deltaT/checkpoint;
    dot = -Pdiss(ite);
        
    plot(ttt(k),dot,'.','MarkerSize',50)
end
hold off

%plot difference between power from viscous disspation and surface deformation
if plotDIFF==1
figure
plot(time,Ps+Pdiss)
ylabel('\Delta E')
xlabel('time')
grid on

end



