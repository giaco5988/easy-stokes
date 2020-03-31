%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

%physical variables
Ca = '6.000000';
lambda = '1.000000';

%options
plotShape = 1; plotMesh = 1;
plotCurv = 0;
plotNormalVectors = 0;

%ID
%D = [0.032 0.0321 0.0322];
D = 0.01;
elem = '500';
dt = 0.001;
totLoop = '3000';
RK = '2';
CPUs = '1';

for  k = 1:numel(D)

    dt = num2str(dt);

    %path
    here = pwd;
    path = '~/Documents/C++/rising_droplet/results/';
    if numel(num2str(D(k)))==6
        results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' num2str(D(k)) '00_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    elseif numel(num2str(D(k)))==5
        results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' num2str(D(k)) '000_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    elseif numel(num2str(D(k)))==4
        if numel(num2str(dt))==5
            results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' num2str(D(k)) '0000_elem=' elem '_dt=' dt '000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
        elseif numel(num2str(dt))==4
            results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' num2str(D(k)) '0000_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
        end
    elseif numel(num2str(D(k)))==3
        results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' num2str(D(k)) '00000_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    else
        warning('No results are loaded');
        break;
    end
    step = 1;

    dt = str2double(dt);
    %dt =0.01;

    cd(results)

    %load options
    A = importdata('ParametersOptions.txt');
    checkpoint = A.data(9);

    loop = str2double(totLoop)/checkpoint+1;

    %allocate memory
    Area = zeros(loop,1);
    V = zeros(loop,1);
    elong = zeros(loop,1);

    for i=1:step:loop

        display([num2str(i) ' of ' num2str(loop)])

        name = ['drop' num2str(i-1) '.txt'];

        try
            A = importdata(name);
            YesBreak = 0;
        catch
           warning('No data')
           YesBreak = 1;
           break;
        end

        x = A.data(:,1); y = A.data(:,2);

        %compute center of mass
        xcm = center_mass(x,y);

        %compute and save area
        Area(i) = surf_gauss_vect(x',y');

        %compute and save volume
        V(i) = axis_int_gauss_vect(x',y');

        %compute rear elongation (positive when tail, negative when indentation)
        elong(i) = x(1)-xcm;

        if plotShape==1
        %plot interface
        figure(1)
        plot([y; -flip(y)],-[x; flip(x)]+xcm,'-')
        if plotMesh==1
            hold on
            plot([y; -flip(y)],-[x; flip(x)]+xcm,'o')
            hold off
        end
        axis equal
        grid on
        xlabel('x')
        ylabel('y')
        cut = 1.5;
        axis([-cut cut -cut cut])
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

    cd(here)

    if YesBreak==1

       V = V(1:i-1);  Area = Area(1:i-1); elong = elong(1:i-1); 

    end

    %plot area variation norm and volume error
    R = nthroot(V(1)/4/pi*3,3);
    A0 = 4*pi*R^2;
    myNorm = sqrt((Area-A0)/(Area(1)-A0));
    figure(5)
    subplot(2,1,1)
    plot(0:dt*checkpoint:dt*(numel(myNorm)-1)*checkpoint,myNorm)
    grid on
    ylabel('||A||')
    xlabel('t')
    hold on
    drawnow

    %plot colume error
    errV = (V-V(1))/V(1);
    subplot(2,1,2)
    plot(0:dt*checkpoint:dt*(numel(myNorm)-1)*checkpoint,errV)
    grid on
    ylabel('err_V')
    xlabel('t')
    hold on
    drawnow

    %plot rear drop elongation
    figure(6)
    plot(0:dt*checkpoint:dt*(numel(myNorm)-1)*checkpoint,(elong-R)/R)
    grid on
    ylabel('rear elongation')
    xlabel('t')
    hold on
    drawnow

end











