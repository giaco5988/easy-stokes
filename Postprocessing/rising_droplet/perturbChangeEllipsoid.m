%check change in pertubation for ellipsoidal droplet

close all
clear variables

%parameters
D = 0.2;
n = 1000;
theta = 0:pi/n:pi;

[x,y]=draw_circle_lean(0,0,n,D);

%plot drop
figure
plot(x,y)
grid on
axis equal
xlabel('x')
ylabel('y')

R = sqrt(x.^2+y.^2);
R0 = 1;
%V0 = axis_int_gauss_vect(R.*cos(theta),R.*sin(theta));
V0 = axis_int_gauss_vect(x,y);
%Van = 4/3*pi;

%errAN = (V0-Van)/Van

figure
plot(theta,R)
xlabel('theta')
ylabel('R')
grid on

%modify perturbation amplitude and check volume conservation
delta = [1e-5 1e-4 1e-3 1e-2 1e-1];
V = zeros(1,numel(delta));
%Rnew = R;
for i = 1:numel(delta)
    
    Rnew = R;
    display([num2str(i) ' of ' num2str(numel(delta))])
    
    Rnew = (Rnew-R0)/max(Rnew-R0)*delta(i) + 1;
    
    figure(3)
    hold on
    plot(Rnew.*cos(theta),Rnew.*sin(theta))
    grid on
    xlabel('x')
    ylabel('y')
    axis equal
    
    figure(2)
    hold on
    plot(theta,Rnew)
    grid on
    xlabel('theta')
    ylabel('R')
    
    V(i) = axis_int_gauss_vect(Rnew.*cos(theta),Rnew.*sin(theta));
    
end

%plot volume error
figure
errV = (V-V0)/V0;
loglog(delta,errV,'o-')
xlabel('\delta')
ylabel('errV')
grid on
title('error on volume')

%compute the analytical error at first order
a = (1+D)^(2/3)/(1-D)^(2/3);    b = (1-D)^(1/3)/(1+D)^(1/3);    %major and minor axis
A = a^2-b^2;    B = b^2;    C = sqrt(B/A);  %dummy variables
errV1 = 3/4*pi*sqrt(A)*(sqrt(1+C^2)+C^2*log(abs(1+sqrt(1+C^2))) + ...
    sqrt(1+C^2)-C^2*log(abs(-1+sqrt(1+C^2)))) -1/3*pi;
errV1 = delta*errV1;

hold on
loglog(delta,errV1,'o-')
hold off
legend('Numerical','Analytical 1st Order','Location','Best')











