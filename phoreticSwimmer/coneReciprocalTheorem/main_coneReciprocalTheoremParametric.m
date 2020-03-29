%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%parameter
manyTheta = linspace(0,pi,101);

%initialize
Unum = zeros(numel(manyTheta),1);
Urec = zeros(numel(manyTheta),1);
for i = 1:numel(manyTheta)

%physical parameters
r = 1;                      % average cone radius
L = 1.8;                   % motor lenght
h = 0.2;                    % wall thickness
theta = manyTheta(i);               % inclinitation of the cone
PARAM.phi = [0 0 0 0];      % concentration
PARAM.flux = [0 0 1 0];     % BC for flux

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;

%geometry parameters
nPerLenght = 50;
PARAM.panels = 4;
PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght)];                  % number of element per panel
%PARAM.typePanel = [0 1 0 1];                % 0 is a straight line, 1 ia an arc

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
xC = L/2;   yC = r;
[x1,y1] = buildStraightLine(L,r+h/2,0,r+h/2,PARAM,1);
[x2,y2] = buildArc(h/2,0,r,[pi/2 3*pi/2],PARAM,2);
[x3,y3] = buildStraightLine(0,r-h/2,L,r-h/2,PARAM,3);
[x4,y4] = buildArc(h/2,L,r,[-pi/2 pi/2],PARAM,4);
x{1} = x1;  x{2} = x2;  x{3} = x3;  x{4} = x4;
y{1} = y1;  y{2} = y2;  y{3} = y3;  y{4} = y4;

%solid body rotation
[x,y] = rigidBodyRotation(x,y,xC,yC,theta,4);

%filename
PARAM.filename = ['coneReciprocalTheorem_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

%upload velocity of moving cone
upName = [PARAM.res '/conePhoreticAllMob_f1=' num2str(PARAM.flux(1)) '_f2=' num2str(PARAM.flux(2)) '_f3=' num2str(PARAM.flux(3)) '_f4=' num2str(PARAM.flux(4)) '_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
%upName = [PARAM.res '/conePhoretic_n=' num2str(sum(PARAM.n)) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
upload = load(upName);
U0 = upload.yStokes(end);
vSlip1 = upload.vSlip1;
vSlip2 = upload.vSlip2;
vSlip3 = upload.vSlip3;
vSlip4 = upload.vSlip4;
% vSlip1 = zeros(PARAM.n(1),1);
% vSlip2 = zeros(PARAM.n(2),1);
% vSlip3 = upload.vSlip;
% vSlip4 = zeros(PARAM.n(4),1);
% l1 = upload.l1;
% l2 = upload.l2;
% l3 = upload.l3;
% l4 = upload.l4;
%l = [l1; l2+l1(end); l3+l1(end)+l2(end); l4+l2(end)+l3(end)+l1(end)];
clear upload

%underlying flow
PARAM.Uunder = -U0;
Unum(i) = U0;

%print to screen
if i==1
    printToScreenStokes(PARAM)
    %printToScreenLaplace(PARAM)
end

%solve Stokes equation
[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

%compute force on the wall
fx = yStokes(1:2:end-1);
fy = yStokes(2:2:end);
F = forceOnPanel(x,y,fx,1:4,PARAM);

%velocity from the reciprocal theorem
[Vrec,Frec] = bodyVelocityReciprocalTheorem(x,y,F,U0,[vSlip1; vSlip2; vSlip3; vSlip4],fx,fy,1:4,PARAM);
%[VrecFlux,FrecFlux] = bodyVelocityReciprocalTheoremFlux(x,y,F,U0,fx,fy,1:4,PARAM);

display(['F from integration is ' num2str(F)])
display(['V from reciprocal theorem is ' num2str(Vrec)])
display(['F from reciprocal theorem is ' num2str(Frec)])
%display(['V from reciprocal theorem FLUX is ' num2str(VrecFlux)])
%display(['F from reciprocal theorem FLUX is ' num2str(FrecFlux)])

Urec(i) = Vrec;

% cd(PARAM.res)
% save(PARAM.filename)
% cd(PARAM.here)

end

%compare velocity from the numerics and from the recirpocal theorem
figure
plot(manyTheta,Unum)
hold on
plot(manyTheta,Urec)
xlabel('\theta')
ylabel('U_{motor}')
title('Motor velocity')
legend('From full numerics','From recirpocal theorem','Location','Best')
grid on






