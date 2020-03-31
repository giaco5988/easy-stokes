%post processing of micromotor

%tic

close all
% clear variables
% 
% %physical variables
% visc = 1;
% Ca = 0.001;
% alpha = 1.6;
% start = 4;
% elem = 570;
% dt = 0.0001;
% theta = -0.05;
% mass = 1;
% 
% if mass==0
%     Ca=1;
% end
% 
% here = pwd;
% res = [here '/LabFrame_mass=' num2str(mass) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=10_alpha=' num2str(alpha) '_InPos=' num2str(start) '_RK2.mat'];
% load(res)

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);
addpath('~/Documents/MATLAB/test/2D_axisymmetric_zeroForce')

m = PARAM.m; q = PARAM.q;
[~,loop] = size(risa);
vel = zeros(1,loop);
mesh = zeros(1,loop);
V = zeros(1,loop);
A = zeros(1,loop);
forceX = zeros(1,loop);
comp = zeros(1,loop);

%stop at this loop
stop = 1700;
step = 10;
plotRealTIME = 1;

%volume at first iteration
indXY = find(risb(:,1)==0,2,'first');
a_drop = risa(m+2:indXY(2)-1,1);
b_drop = risb(m+2:indXY(2)-1,1);
Vin = axis_int_gauss_vect(a_drop',b_drop');

for i = 1:step:loop
    
    disp(i)
    
    indVel = find(risy(:,i)==0,2,'first');
    indXY = find(risb(:,i)==0,2,'first');
    vel(i) = risy(indVel(2)-1,i);
    mesh(i) = indXY(2)-1;
    a_wall = risa(1:m+1,i);
    b_wall = risb(1:m+1,i);
    
%     y = main_noBubble(vel(i));
%     
%     fx = y(1:2:2*PARAM.m);
%     dl = sqrt(diff(b_wall).^2 + diff(a_wall).^2);
%     R = (b_wall(1:end-1)+b_wall(2:end))/2;
%     forceX(i) = 2*pi*sum(fx.*R.*dl);
    
    a_drop = risa(m+2:indXY(2)-1,i);
    b_drop = risb(m+2:indXY(2)-1,i);
    
    phi = FigInOut(a_drop',b_drop',a_wall,b_wall);
    comp(i) = max(phi>pi);
    V(i) = axis_int_gauss_vect(a_drop',b_drop');
    A(i) = surf_gauss_vect(a_drop',b_drop');
    
    %compute excess surface energy
    %DeltaE = 
    %compute power due to dissipation
    
    if plotRealTIME==1
        
%     figure(2)
%     plot(a_wall,b_wall,'b-o',a_wall,-b_wall,'bo-')
%     hold on
%     plot(a_drop,b_drop,'ro-',a_drop,-b_drop,'ro-')
%     hold off
%     axis equal
%     grid on
%     xlabel('x')
%     ylabel('r')
%     axis([-5 11 -2 2])
%     title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc)])
    
    figure(1)
    subplot(3,2,1)
    plot(a_wall,b_wall,'b-o',a_wall,-b_wall,'bo-')
    hold on
    plot(a_drop,b_drop,'ro-',a_drop,-b_drop,'ro-')
    hold off
    axis equal
    grid on
    xlabel('x')
    ylabel('r')
    axis([-5 11 -2.5 2.5])
    title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc)])
    
    %ind = find(risa(:,i)==0, 2, 'first');
    %q = ind(2)-3-m;
%     
    
    
    %plot motor velocity
    %figure(2)
    subplot(3,2,2)
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,vel(i),'ob')
    drawnow
    grid on
    hold off
    xlabel('time')
    ylabel('velocity')
    
    %check compenetration
    subplot(3,2,3)
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,comp(i),'ob')
    grid on
    hold off
    xlabel('time')
    ylabel('InOut')
    
    %plot droplet curvature
    subplot(3,2,4)
    %hold on
    plot(a_drop,KKK(1:indXY(2)-1-m-1,i),'-bo')
    grid on
    %hold off
    xlabel('x')
    ylabel('K')
    
    subplot(3,2,5)
    %hold on
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,A(i),'ob')
    grid on
    hold off
    %hold off
    xlabel('x')
    ylabel('A')
    
    if PARAM.continuity==0
        subplot(3,2,6)
        hold on
        errV = (V(i)-Vin)/Vin;
        plot(i*PARAM.checkpoint*PARAM.deltaT,errV,'ob')
        grid on
        hold off
        xlabel('x')
        ylabel('errV')
    else
        subplot(3,2,6)
        hold on
        plot(i*PARAM.checkpoint*PARAM.deltaT,V(i),'ob')
        grid on
        hold off
        xlabel('x')
        ylabel('V')
    end
    
    end
    
    if i==stop
        disp('STOP ARTIFICIALLY!')
        break;
    end

    
end

%plot first and final configuration
a1wall = risa(1:PARAM.m+1,1);   b1wall = risb(1:PARAM.m+1,1);
a1drop = risa(PARAM.m+2:PARAM.m+PARAM.q+2,1);   b1drop = risb(PARAM.m+2:PARAM.m+PARAM.q+2,1);
a2wall = risa(1:PARAM.m+1,stop);   b2wall = risb(1:PARAM.m+1,stop);
a2drop = risa(PARAM.m+2:PARAM.m+PARAM.q+2,stop);   b2drop = risb(PARAM.m+2:PARAM.m+PARAM.q+2,stop);
figure
subplot(3,1,1)
plot(a1wall,b1wall,'-b',a1wall,-b1wall,'b-','LineWidth',2)
hold on
plot(a1drop,b1drop,'-b',a1drop,-b1drop,'b-','LineWidth',2)
plot(a2wall,b2wall,'-r',a2wall,-b2wall,'r-','LineWidth',2)
plot(a2drop,b2drop,'-r',a2drop,-b2drop,'r-','LineWidth',2)
axis equal
grid on
xlabel('x')
ylabel('r')
title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc)])
hold off

%plot motor velocity
subplot(3,1,2)
hold on
plot((1:step:stop)*PARAM.checkpoint*PARAM.deltaT,vel(1:step:stop),'LineWidth',2)
grid on
xlabel('time')
ylabel('velocity')
hold off

%plot compenetration
subplot(3,1,3)
hold on
plot((1:step:stop)*PARAM.checkpoint*PARAM.deltaT,comp(1:step:stop),'LineWidth',2)
grid on
xlabel('time')
ylabel('InOut')
hold off

%plot mesh elements
% figure
% hold on
% plot((1:stop)*PARAM.checkpoint*PARAM.deltaT,mesh(1:stop),'or')
% grid on
% xlabel('time')
% ylabel('num elem')
% hold off

%compute displacemenet from real position of motor
%(only good for absolute frame of reference)
initPOS = PARAM.L*cos(PARAM.theta);
finalPOS = max(a_wall);
DELTAl = finalPOS-initPOS;
display(DELTAl)

%compute energy of droplet and viscous dissipation
%surface tension based on mass injection NO SENSE!!!!!
radius = PARAM.alpha*PARAM.R; %droplet radius
U = PARAM.mass/4/pi/radius^2;
gamma = abs(U)*PARAM.visc2/PARAM.Ca;

%surface energy
E2 = A(stop)*gamma;
E1 = A(1)*gamma;
deltaE = E2-E1;

%viscous dissipation FROM NUMERICAL SOLUTION
Vavg = sum(vel)/loop;
%DT = PARAM.deltaT*loop;
%Favg = sum(forceX)/loop;
DRAGnum = sum(vel*PARAM.deltaT.*forceX);

%viscous dissipation from analtyical for a pipe
numer = 2*pi*PARAM.visc2*vel;
denom = log(PARAM.L/PARAM.R) - 0.5;
DRAGan = sum(PARAM.deltaT*vel*step*PARAM.checkpoint.*numer/denom);

%display(deltaE)
%display(DRAGnum)
%display(DRAGan)

%displaycement drom numeric
%lNUM = sum(vel*PARAM.deltaT*PARAM.checkpoint*step);

%displacement in Li
xcm = center_mass(a_drop,b_drop);
Rb = sqrt((a_drop(1)-xcm)^2);
L = PARAM.L;
Rj = PARAM.R;
lPAP = 6*Rb^2/(3*Rb+L/(log(2*L/Rj)-0.72));

%anatyical drag used on the paper
Fbubble = -6*pi*Rb*Vavg;
Frod = numer/denom;

%toc

