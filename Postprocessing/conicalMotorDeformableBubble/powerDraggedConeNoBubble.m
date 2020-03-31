% compute velocity of the motor and bubble once the geomtry is known

function [P,F] = powerDraggedConeNoBubble(U,tParametric,PARAM)

PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

disp(['Compute power needed to move an empty cone with velocity U=' num2str(U)])

%modify parameters
PARAM.n = [numel(tParametric{1})-1 numel(tParametric{2})-1 numel(tParametric{3})-1 numel(tParametric{4})-1];
PARAM.panels = 4;
PARAM.blockType = 0;
PARAM.panelType = [0 0 0 0];
PARAM.typeBCstokes = [0 0 0 0];
PARAM.velBCaxial= {0 0 0 0};
PARAM.velBCradial= {0 0 0 0};

%build shape
[xBlock1,yBlock1] = buildGeometryPanelsParametric(tParametric,PARAM);
x = xBlock1(1:4);
y = yBlock1(1:4);

%add underlying flow
PARAM.addFlow = 1;
PARAM.Uunder = U;

%solve Stokes equation
printToScreenStokes(PARAM)
yStokes = BEM_Stokes(x,y,PARAM);

%compute force acting on cone
[xHere,yHere] = getBlockCoordinates(x,y,PARAM,1);
weight = integrationOnLineWeightAxis(xHere,yHere,PARAM.orderVariable(1),PARAM.orderGeometry(1),PARAM.SPlinesType(1));
fx = yStokes(1:2:end-1);
F = weight*fx;

%compute power
P = F*U;

%compute velocity field
xGrid = linspace(-5,15,50);
yGrid = linspace(0,5,50);
[Xgrid,Ygrid] = meshgrid(xGrid,yGrid);
[X0,Y0,uxGrid,uyGrid] = computeVelPressField(Xgrid,Ygrid,x,y,yStokes,0,U,PARAM,0,1,1);

figure
plotFieldAxis(x,y,X0,Y0,uxGrid,uyGrid,1.5,0.2,1,0,PARAM)
axis equal
drawnow
title('Dragged cone')










