%predict motor motion with model

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

close all
clear variables

%problem parameters
%source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/parametric_study/theta005/';
source = '~/Documents/MATLAB/droplet_simulations/server/';
dest = '~/Documents/research_notes/someMovies/motor/frames/';
theta = -0.02;
Ca = 0.001;
dt = 0.01;
visc = 0.1;
inflate = 0;
element = 228;
Inpos = 4;

%options
FrameLab = 1;               % 0 plot in motor frame, 1 lab frame or 2 bubble frame
meshFine = 100;

%upload data
filename = [source 'ConicalMotor_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(element) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=' num2str(Inpos) '_RK=2.mat'];
load(filename)

%solution at this iteration
ite = 2000;

%often used
m = find(risa(2:end,ite)==risa(1,ite));
indNode = find(risb(:,ite)==0,2,'first');
q = indNode(2)-m-3;

%motor and droplet coordinates
aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);

%select solution
numSing = numel(aDrop)+numel(aMotor)-1;
solution = risy(1:2*numSing+1,ite);
velDropX = solution(2*m+1:2:end-2);
velDropY = solution(2*m+2:2:end-1);

%spline coefficients
[ax,bx,cx,dx,ay,by,cy,dy] = spline_symmetric(aDrop',bDrop');

%compute velocity field in the motor oulet
%x = max(aMotor)-0.1;
x = 4;
[~,ind] = min(abs(aMotor(1:round(m/2))-x));
yMax = bMotor(ind);
yMin = 0;
y = linspace(yMin,yMax,meshFine*(yMax-yMin));
[X,Y] = meshgrid(x,y);
[U,V,uFrame] = velocityFieldMicromotor(aDrop',bDrop',aMotor',bMotor',solution,visc,X,Y,PARAM,FrameLab);
U = U-uFrame;

%plot location where velocity is computed and velocity profile
figure
subplot(2,1,1)
plot([aDrop; flip(aDrop)],[bDrop; -flip(bDrop)],'r')
hold on
plot(aMotor,bMotor,'k')
plot(aMotor,-bMotor,'k')
plot(X,Y,'m')
plot(X,-Y,'m')
xlabel('x')
ylabel('r')
axis equal
grid on
axis([-2 12 -2 2])

subplot(2,1,2)
plot(U,Y)
xlabel('U')
ylabel('r')
grid on

%PART FROM SHEAR
%compute average velocity in the outlet
U0 = U(1);    Uavg = U0/2;

%compute film thickness, with the assumptions that it is constant
dist = distWallDrop(aDrop',bDrop',aMotor',bMotor');
h = min(dist);
look = dist<(5*h);
ind1 = find(look==1,1,'first');
ind2 = find(look==1,1,'last');
xFilm = aDrop(ind1:ind2);   yFilm = bDrop(ind1:ind2);

%compute velocity tangent to the interface
%compute normal for nodes
N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2);...
   -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end)).^2+(by(end)+2*cy(end)+3*dy(end)).^2)];
uNormX = velDropX'.*N(1,:);  uNormY = velDropY'.*N(2,:);
uTangX = velDropX'-uNormX;   uTangY = velDropY'-uNormY;
uTangXFilm = uTangX(ind1:ind2);     uTangYFilm = uTangY(ind1:ind2);

%check how is the velocity in the film
velFimlX = velDropX(ind1:ind2);
velFimlY = velDropY(ind1:ind2);
figure
subplot(2,1,1)
%plot(xFilm,yFilm,'r')
plot([aDrop; flip(aDrop)],[bDrop; -flip(bDrop)],'k')
hold on
plot(aMotor,bMotor,'k')
plot(aMotor,-bMotor,'k')
quiver(xFilm,yFilm,velFimlX,velFimlY,'r')
quiver(xFilm,-yFilm,velFimlX,-velFimlY,'r')
xlabel('x')
ylabel('r')
axis equal
grid on

subplot(2,1,2)
%plot(xFilm,velFimlX)
absValue = sqrt(uTangXFilm.^2+uTangYFilm.^2);
Us = absValue.*(uTangXFilm>0)-absValue.*(uTangXFilm<0);
plot(xFilm,Us)
grid on
xlabel('x')
ylabel('U_s')

%mass balance in order to find Us, the velocoty of tanhent to the interface
%Us = Uavg*pi*yMax^2/h;

%shear stress in the film thickness
Shear = Us./dist(ind1:ind2);

%integrate shear along the film lenght
dl = sqrt(diff(xFilm).^2+diff(yFilm).^2);
forceFilm = 2*pi*trapz([0; cumsum(sqrt(diff(xFilm).^2+diff(yFilm).^2))],Shear'.*yFilm);
forceFilmX = forceFilm*cos(theta)

%PART FROM PRESSURE
%compute curvature
K = curv_spline2(bx,by,cx,cy,dx,dy);
Kfront = K(1);  Krear = K(end);

%compute pressure difference between rear and front
deltaPfront = 2*Kfront;
deltaPrear = 2*Krear;
deltaP = 2*(Krear-Kfront);

%integrate the pressure on the circle, rear an front
yFront = bDrop(1:ind1);     yRear = bDrop(ind2:end);
ForceFront = 2*pi*trapz(yFront,yFront)*deltaPfront;
ForceRear = -2*pi*trapz(yRear,yRear)*deltaPrear;

ForceDiff = ForceRear-ForceFront











