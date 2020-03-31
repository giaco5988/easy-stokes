%make graph for transient paper

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

delta = 0.0065;
G = 0.5;

%times
%t = [0  1.22 2.19 3.16 4.13 5.1 6.07 7.04 8.01];
%t = [0  0.55 1.49 2.43 3.37 4.31 5.25 7.13 9.01];
%t = [0 0.825 1.65 2.475 3.3 4.95 6.6 7.425 8.25];
%t = [0 1.3375 2.675 4.0125 5.35 6.6875 8.025 9.3625];
%t = [0 1 10 20 24];
%t = [0 10 14.4 20 30];
t = [0 10 20 23.4 35];
lin_area = zeros(1,numel(t));
lin_vol = zeros(1,numel(t));
nl_vol = zeros(1,numel(t));

%first position of the center of mass;
shift = 9.5;
vert = 2.5;

%figure
hold on

for i=1:numel(t)
    
    disp(i)
    
%     filename = ['crd_rzlam0.5_ca6_mA1000_delta' num2str(delta) '_amp' num2str(G) '_time' num2str(t(i)) '.dat'];
%     A = importdata(filename);
%     x = A(:,1);
%     y = A(:,2);
    
    ite = round(t(i)/deltaT);
    if ite==0; ite=1; end
    xcm = center_mass(risa(:,ite),risb(:,ite));
    plot(risb(1:nbrel(ite)+1,ite)+shift,-risa(1:nbrel(ite)+1,ite)+xcm+vert,'-')
    plot(-risb(1:nbrel(ite)+1,ite)+shift,-risa(1:nbrel(ite)+1,ite)+xcm+vert,'b-')  
    %plot(-x+shift,-y+vert,'r--')
    axis equal
    
    vert = vert+2.5;
    
%     lin_area(i) = surf_gauss_vect(y',x');
%     lin_vol(i) = axis_int_gauss_vect(y',x');
%     nl_vol(i) = V(ite)*V_in;
    
end

hold off

% figure
% plot(deltaT*(0:1:end_loop),Area(1:end_loop+1),'-',t,lin_area,'ro-')
% grid on
% xlabel('t')
% ylabel('Area')
% 
% resize = (nl_vol./lin_vol).^(2/3);
% 
% figure
% plot(deltaT*(0:1:end_loop),Area(1:end_loop+1),'-',t,lin_area.*resize,'ro-')
% grid on
% xlabel('t')
% ylabel('Area resized')