%plot comparison between linear and non linear simulations at different
%time step

%clear all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

%time = [0 0.1 0.2  0.4 1.95 3.5 5.05 6.6 8.15 9.7];
%time = [0 0.2 1.44 3.92 7.64 10.12];
%time = 1:10:loop;
time = 1:1:end_loop;
%time = [0 0.525 1.05 1.575 2.1 3.3 4.5 5.7 6.9 8.1 9.3];
%time = [0 1.325 2.65 3.975 5.3 6.625 7.95 9.275];
%time = [0];
%time = [0 0.825 1.65 2.475 3.3 4.125 5.775 6.6 7.425 8.25 9.075 9.9];
%time = [0 0.825 1.65 2.475 3.3 4.125 5.775 6.6 7.425];


elon_mine = zeros(numel(time),1);
myarea = zeros(numel(time),1);
V_in = axis_int(risa(:,1)',risb(:,1)');

%figure

for i = time
    
    elon_mine(i) = max(risa(1:nbrel(i),i))-min(risa(1:nbrel(i),i));
    %elon_lailai(i) = max(y)-min(y);
    
    %xcm = center_mass(risa(1:nbrel(i)+1,i)',risb(1:nbrel(i)+1,i)');
    
    figure(1)
    plot(risa(1:nbrel(i)+1,i),risb(1:nbrel(i)+1,i),'-k',risa(1:nbrel(i)+1,i),-risb(1:nbrel(i)+1,i),'-k','LineWidth',2)
    axis([-max(max(risa(:,1)),max((risb(:,1))))-0.1 max(max(risa(:,1)),max((risb(:,1))))+0.1 -max(max(risa(:,1)),max((risb(:,1))))-0.1 max(max(risa(:,1)),max((risb(:,1))))+0.1])
    axis equal
    %title(strcat('Comparison between linear and non linear simulation for t=',num2str(time(i)),...
         %' \lambda=',num2str(visc),' Ca=',num2str(Ca)))
    title(strcat('t=',num2str(i*deltaT)))
    hold on
    %plot(-x,-y,'--r',x,-y,'--r','LineWidth',2)
    hold off
    
    myarea(i) = Area(loop);
    
end

% b = a_lailai-4*pi;
% bmod = a_lailai.*resizeV-4*pi;
% myb = myarea-4*pi;
% 
% d = elon_lailai-2;
% myd = elon_mine-2;

% elongation = zeros(round(time(end)/deltaT),1);
% 
% for i = 1:time(end)/deltaT
%     elongation(i) = max(risa(1:nbrel(i),i))-min(risa(1:nbrel(i),i));
% end
% 
% elongation = elongation-2;
% el_norm = sqrt(elongation/d(1));

% figure
% semilogy((1:round(time(end)/deltaT))*deltaT,el_norm,'LineWidth',2)
% legend('mine')
% xlabel('time')
% ylabel('my elongation norm')
% 
% figure
% semilogy(time,myd/myd(1),'o-b',time,d/d(1),'o-r','LineWidth',2)
% legend('mine','Lailai')
% xlabel('time')
% ylabel('elongation norm')
% 
% figure
% plot((0:loop-1)*deltaT,Area(1:loop),'-b',time,a_lailai,'.-r',time,4*pi*ones(numel(time),1),'k-','LineWidth',2,'MarkerSize',25)
% legend('Non linear','Linear','4\pi','LOcation','Best')
% xlabel('time')
% ylabel('Area')
% 
% figure
% plot(time,sqrt(myb/myb(1)),'.b-',time,sqrt(b/b(1)),'.r-','LineWidth',2,'MarkerSize',25)
% legend('Non linear','Linear','LOcation','Best')
% xlabel('time')
% ylabel('A norm')
% 
% figure
% plot((0:loop-1)*deltaT,Area(1:loop),'-b',time,a_lailai.*resizeV,'.-r',time,4*pi*ones(numel(time),1),'k-','LineWidth',2,'MarkerSize',25)
% legend('Non linear','Linear modified','4\pi','LOcation','SouthWest')
% xlabel('time [s]')
% ylabel('Area resized')
% 
% figure
% plot((1:loop-1)*deltaT,V(1:loop-1),'-b',time,v_lailai/v_lailai(1),'.-r','LineWidth',2,'MarkerSize',25)
% legend('Non linear','Linear','LOcation','Best')
% xlabel('time')
% ylabel('Volume')
% 
% figure
% plot(time,sqrt(myb/myb(1)),'.b-',time,sqrt(bmod/bmod(1)),'.r-','LineWidth',2,'MarkerSize',25)
% legend('Non linear','Linear modified','LOcation','Best')
% xlabel('time')
% ylabel('A norm resized')