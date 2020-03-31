%upload edge shape for rising droplet and plot velocity and pressure field

clear variables;
close all;

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%paths
here = pwd;
results = '~/Documents/MATLAB/droplet_simulations/results';
pathUpload = '~/Documents/C++/edge_state/resultsCa=0.25/';

%addpath
addpath('~/Documents/MATLAB/mex_files')

BemDir = '~/Documents/MATLAB/droplet_simulations/drop_buoyancy';

%visualization grid parameters
zoom = 0;
drop_frame = 1;

if zoom==1
    %parameters fot tail visu
    xRange = 0.02; yRange = 0.02;
    shift = 1.01;
    MeshFine = 1000;    
    DistReplace = Inf;
else
    %parameters fot drop visu
    xRange = 2.5; yRange = 2.5;
    shift = 0;
    MeshFine = 30;
    DistReplace = Inf;
end
%what to plot
PlotPressure = 1;

% parameters
n = 200;                                % number of elements
Ca = 0.25;                               % capillary number
lambda = 1;                             % viscosity ratio
UploadShape = 1; elemUpload = 200;      % id upload shape from previous results
RemeshUpload = 0;                       % remesh uploaded shape
IDshape = 8;                           % shape to upload
plotShape = 1;                          % if plotting shape

%directory where to upload the results
UploadResults = [pathUpload 'Ca=' num2str(Ca) '0000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=16_Round=3_elem=' num2str(elemUpload) '_dt=0.000500_loop=50000_RK=2_CPUs=1'];

%name of the saving file
filename = ['RisingDroplet_elem=' num2str(n) '_Ca=' num2str(Ca) '_lambda=' num2str(lambda) '.mat'];

%option numerics
PARAM.way_curv = 1;
PARAM.elem = 1;
PARAM.cfunction = 1;
PARAM.lealpoz = 1;

%upload shape
name = [UploadResults '/drop0.txt'];
A = importdata(name);
    
aIN = A.data(:,1)'; bIN = A.data(:,2)';
a = aIN;    b = bIN;
    
%remesh shape
if RemeshUpload==1
        
    distr = choose_distribution(a,b,1,0.01,n,1,7,1);
    [a, b] = remesh_distribution(a,b,distr);
        
end

%fix first and last point
b([1 end]) = [0 eps];

%compute solution\
cd(BemDir)
[y,~,~,~,N] = bem(a,b,n+1,lambda,1,Ca,PARAM.way_curv,PARAM.lealpoz,PARAM.elem,PARAM);
cd(here)
ux = y(1:2:end-1);  uy= y(2:2:end);
Unormal = ux'.*N(1,:) + uy'.*N(2,:);

%compute drop velocity
Vdrop = DropVelocityAxis(a,b,Unormal');

%compute residuals
res = (ux-Vdrop)'.*N(1,:) + uy'.*N(2,:);
res = max(abs(res));
display(res)

%compute fields
[U,V,p,X,Y,xxx,yyy] = velocity_field_buoyancy(a,b,y,lambda,Ca,-xRange+shift,xRange+shift,yRange,MeshFine,MeshFine,DistReplace);

%plot velocity field
figure
if drop_frame==1
    %quiver(Y,-X,V,-U+Vdrop,'b')
    x = X(1,:); y = Y(:,1)';
    [YYY,XXX] = meshgrid(y,x);
    streamslice(YYY,-XXX,V',-U'+Vdrop,'b')
else
    quiver(Y,-X,V,-U,'b')
end
hold on
%contourf(Y,-X,abs(V))
axis equal
axis([-yRange yRange -xRange-shift xRange-shift])
if drop_frame==1
    %quiver(-Y,-X,-V,-U+Vdrop,'b')
    x = X(1,:); y = Y(:,1)';
    [YYY,XXX] = meshgrid(y,x);
    streamslice(-YYY,-XXX,-V',-U'+Vdrop,'b')
else
    quiver(-Y,-X,-V,-U,'b')
end
%contourf(-Y,-X,abs(V))
plot(yyy,-xxx,'k-',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
hold off
ylabel('axial direction')
xlabel('radial direction')
title('Velocity field')

%plot pressure field
if PlotPressure==1
    figure
    %plot pressure inside
    %InOut = FigInOut(a,b,X(:),Y(:));
    %[lin,col] = size(X);
    %InOut = reshape(InOut,lin,col);
    %contourf(Y,-X,p.*(InOut>pi),'b')
    contourf(Y,-X,p,'b')
    hold on
    contourf(-Y,-X,p,'b')
    axis equal
    axis([-yRange yRange -xRange-shift xRange-shift])
    plot(yyy,-xxx,'k-',-yyy,-xxx,'k','MarkerSize',2,'LineWidth',2)
    hold off
    ylabel('axial direction')
    xlabel('radial direction')
    title('Pressure field')
    colorbar
end
















