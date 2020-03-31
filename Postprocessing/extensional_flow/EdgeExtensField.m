%upload edge shape for rising droplet and plot velocity and pressure field

clear variables;
close all;

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%paths
here = pwd;
results = '~/Documents/MATLAB/droplet_simulations/results';
pathUpload = '~/Documents/C++/edge_state_extensional/resultsServer2/';

%addpath
addpath('~/Documents/MATLAB/mex_files')

BemDir = '~/Documents/MATLAB/droplet_simulations/drop_extensional';

%visualization grid parameters
zoom = 0;

if zoom==1
    %parameters fot tail visu
    xRange = 0.02; yRange = 0.02;
    shift = 1.01;
    MeshFine = 1000;    
    DistReplace = Inf;
else
    %parameters fot drop visu
    xRange = 3; yRange = 2;
    shift = 0;
    MeshFine = 5;
    DistReplace = Inf;
end
%what to plot
PlotPressure = 0;

% parameters
n = 200;                                % number of elements
Ca = 0.1;                               % capillary number
lambda = 1;                             % viscosity ratio
UploadShape = 1; elemUpload = 200;      % id upload shape from previous results
RemeshUpload = 0;                       % remesh uploaded shape
IDshape = 5;                            % shape to upload
plotShape = 1;                          % if plotting shape

%directory where to upload the results
UploadResults = [pathUpload 'Ca=' num2str(Ca) '00000_lambda=' num2str(lambda) '.000000_IDshape=' num2str(IDshape) '_IDdelta=16_Round=4_elem=' num2str(elemUpload) '_dt=0.010000_loop=25000_RK=2_CPUs=1'];

%name of the saving file
%filename = ['ExtensionalDroplet_elem=' num2str(n) '_Ca=' num2str(Ca) '_lambda=' num2str(lambda) '.mat'];

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
[y,~,~,~,N] = bem_extens(a,b,n+1,lambda,1,PARAM.way_curv,PARAM.elem,Ca,PARAM);
cd(here)
ux = y(1:2:end-1);  uy= y(2:2:end);
Unormal = ux'.*N(1,:) + uy'.*N(2,:);

%compute drop velocity
Vdrop = DropVelocityAxis(a,b,Unormal');

%compute fields
[U,V,p,X,Y,xxx,yyy] = velocity_field_extensional(a,b,y,lambda,Ca,-xRange+shift,xRange+shift,yRange,MeshFine,MeshFine,DistReplace);

%velocity magnitude
absU = sqrt(U.^2+V.^2);

%plot velocity field
figure
quiver(X,-Y,U,-V,'b')
hold on
axis equal
axis([ -xRange+shift xRange+shift -yRange yRange])
%contourf(X,Y,absU)
streamslice(X,Y,U,V)
plot(xxx,yyy,'k-',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
hold off
ylabel('axial direction')
xlabel('radial direction')
title('Velocity field')
%colorbar

%plot pressure field
if PlotPressure==1
    figure
    contourf(X,Y,p,'b')
    hold on
    contourf(X,-Y,p,'b')
    axis equal
    axis([ -xRange+shift xRange+shift -yRange yRange])
    plot(xxx,yyy,'k-',xxx,-yyy,'k','MarkerSize',2,'LineWidth',2)
    hold off
    ylabel('axial direction')
    xlabel('radial direction')
    title('Pressure field')
    colorbar
end
















