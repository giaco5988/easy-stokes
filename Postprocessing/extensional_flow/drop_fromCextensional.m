%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variables
%manyCa = 0.01:0.01:0.11;
manyCa = 0.05;
lambda = '1.000000';

plotD = 1;

%options
plotShape = 1;
plotCurv = 0;
plotNormalVectors = 0;

%ID
manyD = 0.02;
%manyD = [0.02 0.03 0.05 0.08 0.1 0.13 0.16 0.2 0.24 0.3 0.3];
%manyD = [0.02 0.03 0.08 0.1 0.13 0.16 0.2 0.24 0.3 0.3];
elem = '50';
dt = '0.100000';
totLoop = '500';
RK = '2';
CPUs = '1';

%path
here = pwd;
path = '~/Documents/C++/drop_extensionalMPI/results2/';

%initialize
Dfinal = zeros(1,numel(manyCa));

for l = 1:numel(manyCa)
Ca = manyCa(l);
D = manyD(l);

display(['Ca=' num2str(Ca) ' visc=' lambda])

dt = '0.100000';
    
results = [path 'Ca=' num2str(Ca) '0000_lambda=' lambda '_delta=' num2str(D) '0000_elem=' elem '_dt=' dt '_loop=' totLoop '_RK=' RK ''];
if Ca==0.1
results = [path 'Ca=' num2str(Ca) '00000_lambda=' lambda '_delta=' num2str(D) '00000_elem=' elem '_dt=' dt '_loop=' totLoop '_RK=' RK ''];
end
if D==0.1||D==0.2||(D==0.3&&Ca==0.11)
results = [path 'Ca=' num2str(Ca) '0000_lambda=' lambda '_delta=' num2str(D) '00000_elem=' elem '_dt=' dt '_loop=' totLoop '_RK=' RK ''];
end

dt = str2double(dt);
%dt =0.01;

cd(results)

%load options
A = importdata('ParametersOptions.txt');
checkpoint = A.data(8);

loop = str2double(totLoop)/checkpoint;
step = 1;

%allocate memory
Area = zeros(loop,1);
V = zeros(loop,1);
D = zeros(loop,1);

for i = 1:step:loop
    
    display([num2str(i) ' of ' num2str(loop)])
    
    name = ['drop' num2str(i-1) '.txt'];
    
    A = importdata(name);

    x = A.data(:,1); y = A.data(:,2);
    
    %compute center of mass
    xcm = center_mass(x,y);
    
    %compute and save area
    Area(i) = surf_gauss_vect(x',y');

    %compute and save volume
    V(i) = axis_int_gauss_vect(x',y');
    
    L = max(x);
    B = max(y);
    D(i) = (L-B)/(L+B);
    
    if plotShape==1
    %plot interface
    figure(1)
    plot([x; flip(x)],[y; -flip(y)],'o-')
    axis equal
    grid on
    xlabel('x')
    ylabel('y')
    axis([-2+xcm 2+xcm -2 2])
    drawnow
    end
    
    if plotCurv==1
    figure(2)
    k = A.data(:,3)+A.data(:,4);
    %k = A.data(:,4);
    %plot curvature
    plot(x-xcm,k,'o-')
    grid on
    xlabel('x')
    ylabel('k')
    drawnow
    end
    
    if plotNormalVectors==1
    figure(3)
    subplot(2,1,1)
    nx = A.data(:,5);
    ny = A.data(:,6);
    %plot curvature
    plot(x-xcm,nx,'o-')
    grid on
    xlabel('x')
    ylabel('nx')
    subplot(2,1,2)
    plot(x-xcm,ny,'o-')
    grid on
    xlabel('x')
    ylabel('ny')
    drawnow
    end
    
end

Dfinal(l) = D(i);

if plotD==1
%plot ellipcity and volume error
figure
subplot(2,1,1)
plot(0:dt*step:dt*(numel(D)-1),D(1:step:loop),'k')
grid on
ylabel('D')
xlabel('t')
title(['Ca=' Ca ' \lambda=' lambda])

errV = (V-V(1))/V(1);
subplot(2,1,2)
plot(0:dt*step:dt*(numel(D)-1),errV(1:step:loop),'k')
grid on
ylabel('err_V')
xlabel('t')

%plot D residuals
figure
semilogy(dt:dt*step:dt*(numel(D)-1),diff(D(1:step:loop)))
grid on
ylabel('res_D')
xlabel('t')
title(['Ca=' Ca ' \lambda=' lambda])
end

end

filename = '~/Documents/MATLAB/droplet_simulations/results/extensional/viscosity/';
cd(filename)
file_x2 = strcat('xxx_1.mat');
file_y2 = strcat('yyy_1.mat');
load(file_x2);  load(file_y2);

%plot results
figure
hold on
plot(manyCa,Dfinal,'o-k')
grid on
xlabel('Ca')
ylabel('D')
title(['\lambda=' lambda])
plot(xxx,yyy,'--b','LineWidth',2)

legend('my data c++','Stone et al. 1989','Location','Best')



cd(here)
