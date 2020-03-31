%compare kinetic and surface energy

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%load data
%load('~/Documents/MATLAB/droplet_simulations/results/micromotor/LabFrame_theta=-0.05_el=570_dt=0.01_visc=1_Ca=1_R=1_L=10_alpha=1.6_InPos=5_RK2.mat')
load('~/Documents/MATLAB/droplet_simulations/results/micromotor/LabFrame_theta=-0.05_el=670_dt=0.01_visc=1_Ca=1_R=1_L=10_alpha=1.6_InPos=-2_RK2.mat')

stop = 150;
step = 1;
plotNOW = 0;
plotVEL = 0;

fixed_elem = m;

%select BUBBLE solution and geometry
risaBUB = risa(PARAM.m+2:PARAM.m+PARAM.q+2,:);    risbBUB = risb(PARAM.m+2:PARAM.m+PARAM.q+2,:);    risyBUB = risy(2*PARAM.m+1:2*PARAM.m+3+2*PARAM.q,:);

%geometry at first
aIN = risaBUB(:,1)';   bIN = risbBUB(:,1)';
Rin = PARAM.alpha*PARAM.R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute force on motor alone with U=1, than it scales linearly
% [y,aaa,bbb] = main_noBubble(1,PARAM.theta,PARAM.L,PARAM.R);
% fx = y(1:2:end-1);
% dl = sqrt(diff(bbb).^2 + diff(aaa).^2);
% r = (bbb(1:end-1)+bbb(2:end))/2;
% forceX = 2*pi*sum(fx'.*r.*dl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
Area = zeros(stop+1,1);
Area(end) = surf_gauss_vect(risaBUB(:,stop+1)',risbBUB(:,stop+1)');
xcmBUB = zeros(stop+1,1);
xcmBUB(end) = center_mass(risaBUB(:,stop+1)',risbBUB(:,stop+1)');
magn = zeros(stop,1);
FX = zeros(stop,1);
Prepulse = zeros(stop,1);
visc_diss = zeros(stop,1);
for i = 1:stop
    
    disp([num2str(i) ' of ' num2str(stop)])
    
    %all geometry
    aaa = risa(1:PARAM.m+PARAM.q+2,i)';   bbb = risb(1:PARAM.m+PARAM.q+2,i)';
    
    %bubble geometry
    a = risaBUB(:,i)';
    b = risbBUB(:,i)';
    
    %bubble center of mass
    xcmBUB(i) = center_mass(a,b);
    
    %plot droplet
    if plotNOW==1 && plotVEL==0
        figure(1)
        plot(a,b,'b',a,-b,'b')
        axis equal
        grid on
        %axis([-xlim xlim -ylim ylim])
        xlabel('x')
        ylabel('r')
        drawnow
    end
    
    %compute area
    Area(i) = surf_gauss_vect(a,b);
    
    %plot drolet and velocity field
    if plotVEL==1
        figure(1)
        plot(a,b,'b',a,-b,'b')
        hold on
        quiver(X,Y,U,V,'b')
        quiver(X,-Y,U,-V,'b')
        axis equal
        grid on
        %axis([-xlim xlim -ylim ylim])
        xlabel('x')
        ylabel('r')
        hold off
        drawnow
    end
    
    %compute distance between droplet and wall FOR DROP
    distDWdrop = distDropWall(a,b,aaa(1:PARAM.m+1),bbb(1:PARAM.m+1));
    %ux = risyBUB(1:2:end-1,i)' - risy(2*PARAM.q+2*PARAM.m+3,i);  uy = risyBUB(2:2:end,i)';
    %display(risy(2*PARAM.q+2*PARAM.m+3,i))
    ux = risyBUB(1:2:end-1,i)';  uy = risyBUB(2:2:end,i)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute energy due to repulsive forces
    %add Van der Walls with Hamacker model
%     ForceVanDerWalls = -PARAM.strength*distDWdrop.^(-3);
%     %df_x = df_x + ForceVanDerWalls.*N(1,:);
%     %df_y = df_y + ForceVanDerWalls.*N(2,:);
%     int = ForceVanDerWalls.*N(1,:).*ux + ForceVanDerWalls.*N(2,:).*uy;
%     int1 = int(1:end-1);    int2 = int(2:end);  dl = sqrt(diff(b).^2+diff(a).^2);
%     INT = pi*(int1.*b(1:end-1)+int2.*b(2:end)).*dl;
%     Prepulse(i) = sum(INT);  
%     
%     %compute integral of VanDerWalls in order to balance zero force
%     %condition
%     temp = ForceVanDerWalls.*N(1,:);   dl = sqrt(diff(aaa(fixed_elem+2:2+PARAM.m+PARAM.q)).^2+diff(bbb(fixed_elem+2:2+PARAM.m+PARAM.q)).^2);
%     FX(i) = pi*sum((temp(1:end-1).*bbb(fixed_elem+2:2+PARAM.m+PARAM.q-1)+temp(2:end).*bbb(fixed_elem+3:2+PARAM.m+PARAM.q)).*dl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute viscous dissipation on the surface
    %velocities
    [df_x,df_y] = df_surf_tens_spline_symm(a,b,1/PARAM.Ca);  %BC stresses
    visc_diss(i) = int_axis_spline_symmetric(a,b,df_x'.*ux' + df_y'.*uy');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %compute force on stucK motor
    %y = BEM_Force(aaa,bbb,PARAM);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %compute stresses on cone at rest due to a certain velocity
    %ris = main_noBubble(risy(2+PARAM.m+PARAM.q));
    
%     figure(2)
%     plot(ux,'o-')
%     hold on
%     %plot(uy)
%     %hold off
%     xlabel('x')
%     ylabel('ux')
%     grid on
%     drawnow
    
end

%compute energy
Es = Area-4*pi*Rin^2;   %surface energy
Pdiss = visc_diss/PARAM.Ca;
Ps = -diff(Es)/PARAM.deltaT/PARAM.checkpoint;

time = 0:PARAM.checkpoint*PARAM.deltaT:PARAM.checkpoint*PARAM.deltaT*stop;

figure
plot(time(1:end-1),Ps)
hold on
plot(time(1:step:end-1),-Pdiss(1:step:end),'o')
%plot(time(1:step:end-1),-Pdiss(1:step:end)+risy(2*PARAM.q+2*PARAM.m+3,1:step:stop)'/8/pi,'s-')
%plot(time(1:step:end-1),-Prepulse(1:step:end)-Pdiss(1:step:end)+FX(1:step:stop).*risy(2*PARAM.q+2*PARAM.m+3,1:step:stop)','s-')
hold off
grid on
xlabel('t')
ylabel('Power')
legend('dS/st','Dissipation','Location','Best')

DE = Pdiss+Ps;

figure
subplot(2,1,1)
%plot(time(1:end-1),Ps)
plot(time(1+step:step:end-1),DE(1+step:step:end))
hold on
%plot(time(1+step:step:end-1),risy(2*PARAM.q+2*PARAM.m+3,1+step:step:stop)*forceX)
hold off
grid on
ylabel('\Delta E')
xlabel('t')

subplot(2,1,2)
plot(time(1+step:step:end-1),risy(2*PARAM.q+2*PARAM.m+3,1+step:step:stop))
xlabel('t')
grid on
ylabel('U')
hold off
%legend('\Delta E','U')

%plot center of mass position of the bubble and motor displacement
figure
plot(time,xcmBUB-xcmBUB(1),'o')
hold on
plot(time,-risa(1,1:stop+1)+risa(1,1),'o')
xlabel('time')
ylabel('x_{CM}')
grid on







