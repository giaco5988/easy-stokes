%nonlinear function for motor model

function [R,first,second,third] = nonlinearFunctionModel(x,mu,gamma,chi,V,alpha,Lf,Li,D1,x0)

%compute geometry
%A = findR1R2HforModel(PARAM,Vtot,x,alpha);

%compute derivative
xp = D1*x;

%first term
first = V/alpha./x.^2;

%second
second = 2*pi*chi*(mu/gamma)^(2/3)*nthroot(xp.^2,3);

%third
third = -4*pi*mu/gamma*xp.*x.^2.*(1/Lf-1/Li-4*(V*pi*alpha^2*x.^2)./(pi^2*alpha^4*x.^6-V^2));

%residuals
R = second+third-first;
%R = third-first;
R(1) = x0-x(1);
%R = -third-first;