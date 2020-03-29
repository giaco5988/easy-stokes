%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 0.25;                      % average cone radius
L = 1;                   % motor lenght
h = 0.1;                    % wall thickness
theta = 0;               % inclinitation of the cone
PARAM.fluxBC{1} = 0;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = 0;
PARAM.fluxBC{4} = 0;

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;
compareWithFullSimulation = 0;

%geometry parameters
nPerLenght = 200;
PARAM.panels = 4;
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [0 0 0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 0;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.addFlow = 1;                        % 1 add undelying flow
PARAM.deflationBlock = 1;

%shape function for BCs
PARAM.velBC{1} = 0;
PARAM.velBC{2} = 0;
PARAM.velBC{3} = 0;
PARAM.velBC{4} = 0;

%build geometry
xC = 0;   yC = r;
[x1,y1] = buildStraightLine(L,r+h/2,0,r+h/2,PARAM,1);
[x2,y2] = buildArc(h/2,0,r,[pi/2 3*pi/2],PARAM,2);
[x3,y3] = buildStraightLine(0,r-h/2,L,r-h/2,PARAM,3);
[x4,y4] = buildArc(h/2,L,r,[-pi/2 pi/2],PARAM,4);
x{1} = x1;  x{2} = x2;  x{3} = x3;  x{4} = x4;
y{1} = y1;  y{2} = y2;  y{3} = y3;  y{4} = y4;

%solid body rotation
[x,y] = rigidBodyRotation(x,y,xC,yC,theta,1:4);

%plot conical motor shape
figure
plot(x{1},y{1},'k')
hold on
plot(x{2},y{2},'k')
plot(x{3},y{3},'k')
plot(x{4},y{4},'k')
plot(x{1},-y{1},'k')
plot(x{2},-y{2},'k')
plot(x{3},-y{3},'k')
plot(x{4},-y{4},'k')
grid on
xlabel('l')
ylabel('r')
axis equal
axis([-1.5 1.5 -1.5 1.5])
title(['Motor geometry \theta=' num2str(theta)])

%filename
PARAM.filename = ['coneReciprocalTheorem_n=' num2str(nPerLenght) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

%upload velocity of moving cone
if compareWithFullSimulation==1
    upName = [PARAM.res '/conePhoreticAllMob_f1=' num2str(PARAM.BCflux{1}) '_f2=' num2str(PARAM.BCflux{2}) '_f3=' num2str(PARAM.BCflux{3}) '_f4=' num2str(PARAM.BCflux{4}) '_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
    upload = load(upName);
    U0 = upload.yStokes(end);
    vSlip1 = upload.vSlip1;
    vSlip2 = upload.vSlip2;
    vSlip3 = upload.vSlip3;
    vSlip4 = upload.vSlip4;
    l1 = upload.l1;
    l2 = upload.l2;
    l3 = upload.l3;
    l4 = upload.l4;
    l = [l1; l2+l1(end); l3+l1(end)+l2(end); l4+l2(end)+l3(end)+l1(end)];
    clear upload
elseif compareWithFullSimulation==0
    U0 = 1;
    l = computeArcLenghtBlock(x,y,PARAM,1);
end

%underlying flow
PARAM.Uunder = -U0;

%print to screen
printToScreenStokes(PARAM)

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%plot stresses on wall
figure
fx = yStokes(1:2:end-1);
fy = yStokes(2:2:end);
plot(l,fx)
hold on
plot(l,fy)
grid on
xlabel('x')
ylabel('f_x,f_y')
legend('f_x','f_y','Location','Best')
title('Stress on the boundary')

%plot tangent stresses on wall
figure
%fx = yStokes(1:2:end-1);
%fy = yStokes(2:2:end);
PARAM.orderGeometry = PARAM.orderGeometryStokes;
[~,~,nx,ny,panelRange] = getBlockCoordinates(x,y,PARAM,1);
tx = -ny;   ty = nx;
fTangent = fx.*tx' + fy.*ty';
plot(l,fTangent)
grid on
xlabel('l')
ylabel('t \cdot \sigma \cdot n')
title('Tangent stress from fluid to wall')

%plot tangent stresses on wall, x and y componenet
% figure
% fx = yStokes(1:2:end-1);
% fy = yStokes(2:2:end);
% fTangentX = fTangent.*tx';
% fTangentY = fTangent.*ty';
% plot(l,fTangentX,l,fTangentY)
% grid on
% xlabel('l')
% legend('(t \cdot \sigma \cdot n)_x','(t \cdot \sigma \cdot n)_y','Location','Best')
% ylabel('t \cdot \sigma \cdot n')
% title('Tangent stress')

%compute force on the wall
F = forceOnPanel(x,y,fx,1:4,PARAM);
Ftangent = forceOnPanel(x,y,fTangent,1:4,PARAM);

%output
display(['F from integration is ' num2str(F)])
display(['Integration of tangent force is ' num2str(Ftangent)])

%velocity from the reciprocal theorem
%chabge sign of F because is from fluid to the wall
if compareWithFullSimulation==1
    [Vrec,Frec] = bodyVelocityReciprocalTheorem(x,y,-F,U0,[vSlip1; vSlip2; vSlip3; vSlip4],fx,fy,1:4,PARAM);
    %[VrecFlux,FrecFlux] = bodyVelocityReciprocalTheoremFlux(x,y,F,U0,fx,fy,1:4,PARAM);

    display(['V from reciprocal theorem is ' num2str(Vrec)])
    display(['F from reciprocal theorem is ' num2str(Frec)])
    %display(['V from reciprocal theorem FLUX is ' num2str(VrecFlux)])
    %display(['F from reciprocal theorem FLUX is ' num2str(FrecFlux)])
    
end

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)









