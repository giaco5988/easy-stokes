%try volume correction

close all
clear variables

%draw drop
D = 0.2;
n = 100;
theta = 0:pi/n:pi;

[x,y]=draw_circle_lean(0,0,n,D);

V1 = axis_int_gauss_vect(x,y);
A = surf_gauss_vect(x,y);
%compute the spline coeff
  [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (x, y);
  %compute the versor normal tp the node (N) and to the element (r)
  %[N,r] =normal_versor_clean(a,b,q-1);
N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
      -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

%resize with volume correction
nx = N(1,:);    ny = N(2,:);
Vtarget = 2*V1;
displace = -3*(V1-Vtarget)/A;
zx = displace*nx;
zy = displace*ny;

Xnew = x + zx;  Ynew = y + zy;
V2 = axis_int_gauss_vect(Xnew,Ynew);

figure
plot(x,y)
hold on
axis equal
grid on
plot(Xnew,Ynew)
xlabel('x')
ylabel('y')
legend('original','resized','Location','Best')

%error on target volume
errVtarg = (V2-Vtarget)/Vtarget;
display(errVtarg)