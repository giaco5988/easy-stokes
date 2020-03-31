%model startin with reyssat

clear variables
close all

%time to compute
T = 5e3;
n = 1e2;

%geometrical parameters
alpha = 0.05;
size = 1.2;
V = 4/3*pi*size^3;
L = 10;
radius = 1;
Li = (1-sin(alpha))/tan(alpha);
Lf = Li+L;
x0 = 5+Li;

%hydro parameters
gamma = 1;
chi = 2.47;
mu = 1;

%derivatives
scale = -2/T;
[t,DM] = chebdif(n,2);
t = (t-1)/scale;
D1 = DM(:,:,1)*scale;
D2 = DM(:,:,2)*scale^2;

%nonlinear function for bubble displacement
f = @(x) nonlinearFunctionModel(x,mu,gamma,chi,V,alpha,Lf,Li,D1,x0);

%newton method
options = optimoptions('fsolve','TolFun',1e-10,'TolX',1e-10,'Display','iter','MaxFunEvals',1e4,'MaxIter',400);
x = fsolve(f,Li*ones(numel(t),1),options);
vDrop = D1*x;

%value of the different terms for the solution
[R,first,second,third] = f(x);

FonWall = second+third;
%FonWall = third;
vMotor = FonWall*(log(L/2/radius)+log(2)-0.5)/(2*pi*L);

%plot results
figure
plot(t,first)
hold on
plot(t,second)
plot(t,third)
grid on
xlabel('t')
ylabel('F')
legend('F_{\gamma}','F_{r}+F_{a}','F_b','Location','Best')
title('Forces acting in the CV')

figure
plot(t,x-Li)
grid on
xlabel('t')
ylabel('x')
title('Bubble position')

figure
plot(t,-vDrop)
grid on
xlabel('t')
ylabel('v')
title('Bubble velocity')

% figure
% plot(t,vMotor)
% grid on
% xlabel('t')
% ylabel('v')
% title('Motor velocity')










