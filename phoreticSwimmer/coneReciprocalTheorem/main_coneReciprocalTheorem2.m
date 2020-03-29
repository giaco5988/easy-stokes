%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                      % average cone radius
L = 1.79;                   % motor lenght
h = 0.2;                    % wall thickness
theta = pi/8;               % inclinitation of the cone
PARAM.phi = [0 0 0 0];      % concentration
PARAM.flux = [1 1 1 1];     % BC for flux

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
plotField = 0;

%geometry parameters
nPerLenght = 100;
PARAM.panels = 4;
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel
%PARAM.typePanel = [0 1 0 1];                % 0 is a straight line, 1 ia an arc

%numerics parameters for Stokes
PARAM.typeBCstokes = [4 4 4 4];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity, 4 is prescibed axial velocity
PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [0 0 0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 0;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 0;                        % 0 is zero flow at infinity, 1 add undelying flow

%shape function for BCs
PARAM.velBC{1} = 1;
PARAM.velBC{2} = 1;
PARAM.velBC{3} = 1;
PARAM.velBC{4} = 1;

%build geometry
xC = L/2;   yC = r;
[x1,y1] = buildStraightLine(L,r+h/2,0,r+h/2,PARAM,1);
[x2,y2] = buildArc(h/2,0,r,[pi/2 3*pi/2],PARAM,2);
[x3,y3] = buildStraightLine(0,r-h/2,L,r-h/2,PARAM,3);
[x4,y4] = buildArc(h/2,L,r,[-pi/2 pi/2],PARAM,4);
x{1} = x1;  x{2} = x2;  x{3} = x3;  x{4} = x4;
y{1} = y1;  y{2} = y2;  y{3} = y3;  y{4} = y4;

%solid body rotation
[x,y] = rigidBodyRotation(x,y,xC,yC,theta,4);

%plot conical motor shape
figure
plot(x{1},y{1},'k')
hold on
plot(x{2},y{2},'k')
plot(x{3},y{3},'r')
plot(x{4},y{4},'k')
plot(x{1},-y{1},'k')
plot(x{2},-y{2},'k')
plot(x{3},-y{3},'r')
plot(x{4},-y{4},'k')
grid on
xlabel('x')
ylabel('r')
axis equal
title('Motor geometry')

%filename
PARAM.filename = ['coneReciprocalTheorem_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

%upload velocity of moving cone
% upName = [PARAM.res '/conePhoreticAllMob_f1=' num2str(PARAM.flux(1)) '_f2=' num2str(PARAM.flux(2)) '_f3=' num2str(PARAM.flux(3)) '_f4=' num2str(PARAM.flux(4)) '_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
% upload = load(upName);
% U0 = upload.yStokes(end);
% clear upload

%underlying flow
%PARAM.Uunder = -U0;

%print to screen
printToScreenStokes(PARAM)

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%plot stresses on wall
figure
fx = yStokes(1:2:end-1);
fy = yStokes(2:2:end);
plot(Xsing,fx)
hold on
plot(Xsing,fy)
grid on
xlabel('x')
ylabel('f_x,f_y')

%compute force on the wall
F = forceOnPanel(x,y,fx,1:4,PARAM);

%velocity from he reciprocal theorem
%[Vrec,Frec] = bodyVelocityReciprocalTheorem(x,y,F,U0,fx,fy,1:4,PARAM);

display(['F from integration is ' num2str(F)])
%display(['V from reciprocal theorem is ' num2str(Vrec)])
%display(['F from reciprocal theorem is ' num2str(Frec)])

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









