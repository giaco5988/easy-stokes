%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variables
Ca = '0.250000';
lambda = '1.000000';
oblate = 0;

%options
plotShape = 0;  zoom = 0;   plotMesh = 0;
plotCurv = 0; 
plotNormalVectors = 0;
plotForceFree = 0;
plotVolume = 0;
plotAreaNorm = 0;
step = 10;

colors= 'brgkcm';
colors = repmat(colors,1,1000);

%path
here = pwd;
if oblate==1
    
    ob = '_oblate';
    
else
    
    ob = '';
    
end

%load path
if strcmp(Ca,'1.000000')
    path = ['~/Documents/C++/edge_state/resultsCa=1' ob '/'];
elseif strcmp(Ca,'0.100000')
    path = '~/Documents/C++/edge_state/resultsCa=0.1/';
elseif strcmp(Ca,'0.500000')
    path = '~/Documents/C++/edge_state/resultsCa=0.5/';
elseif strcmp(Ca,'0.250000')
    path = '~/Documents/C++/edge_state/resultsCa=0.25/';
end
%path = '~/Documents/C++/edge_state/resultsBella/';
%path = '~/Documents/C++/edge_state/resultsExtra5/';
%path = '~/Documents/C++/edge_state/resultsOblate2/';
%path = '~/Documents/C++/edge_state/resultsRefineOblate/';
%path = '~/Documents/C++/edge_state/resultsRefineProlate/';
%path = '~/Documents/C++/edge_state/resultsCa=50/';
%path = '~/Documents/C++/edge_state/results4/';
%path = '~/Documents/C++/edge_state/results3/';
%path = '~/Documents/C++/edge_state/results2/';
%path = '~/Documents/C++/edge_state/results500elements/';
%path = '~/Documents/C++/edge_state/resultsDtLess/';
%path = '~/Documents/C++/edge_state/resultsMoreElements/';
%path = '~/Documents/C++/edge_state/results/';

%ID
IDshape = 2:10;    BreakShape = Inf;
IDdelta = 1:16;
Round = 1;
CheckRound = 1;
numIDdelta = numel(IDdelta);
elem = '500';
dt = 0.0005;
totLoop = '50000';
RK = '2';
CPUs = '1';

if numel(IDdelta)==1
    only2 = 1;
else
    only2 = 0;
end

Tstart = 0;

for l = 1:numel(IDshape)
    
    dt = num2str(dt);
    
    if CheckRound==1
        Round = 0;
        findRound = 0;
        %find out rigth round
        while findRound<1000
            
           %dt = str2double(dt);

           findRound = findRound+1;
           if  numel(num2str(dt))==6
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '00_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
           elseif  numel(num2str(dt))==4
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
           elseif numel(num2str(dt))==5
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
           end
           
           cd(results)

           %open .txt
           roundTXT = importdata('round.txt');
           Round = roundTXT(1,1);

           if Round==1
               break;
           end

        end
        
        %findRound = Round;
        
    else
    
        findRound = Round;
    
    end
    %if only2==1
    %find the IDdelta that set the transition
    for m = 1:numIDdelta

        if  numel(num2str(dt))==4
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(m) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
        elseif numel(num2str(dt))==5
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(m) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
        elseif numel(num2str(dt))==6
                results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(m) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '00_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
        end
        cd(results)
        A = importdata('break.txt');
        myBreak = A(1,1);
        if myBreak==1
            break;
        end
        
    end
    %m=1;
    
    %IDdelta = [m-1 m];
    %end
    
    for  k = 1:numel(IDdelta)

    dt = num2str(dt);

    if  numel(num2str(dt))==4
        results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(IDdelta(k)) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    elseif numel(num2str(dt))==5
        results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(IDdelta(k)) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    elseif numel(num2str(dt))==6
        results = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(IDdelta(k)) '_Round=' num2str(findRound) '_elem=' elem '_dt=' dt '00_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
    end

    dt = str2double(dt);
    %dt =0.01;

    cd(results)

    %load options
    A = importdata('ParametersOptions.txt');
    checkpoint = A.data(12);

    loop = str2double(totLoop)/checkpoint+1;

    %allocate memory
    Area = zeros(loop,1);
    V = zeros(loop,1);
    elong = zeros(loop,1);
    deltaNorm = zeros(loop,1);

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
        %display(['xcm=' num2str(xcm)])
        %x(1)

        %compute and save area
        Area(i) = surf_gauss_vect(x',y');

        %compute and save volume
        V(i) = axis_int_gauss_vect(x',y');

        %compute rear elongation (positive when tail, negative when indentation)
        R0 = nthroot(V(i)/4/pi*3,3);
        elong(i) = x(1) - xcm - R0;
        
        %compute norm of perturbation
        deltaNorm(i) = norm(sqrt((x-xcm).^2+y.^2)-R0);

        if plotShape==1
            %compute the spline coeff
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (x', y');

            %compute splines coordinates
            t = 0:0.05:0.9;
            ttt = repmat(t,1,numel(ax));
            axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
            bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
            cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
            dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
            ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
            byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
            cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
            dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));

            %splines coordinates
            xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 x(end)];
            yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 y(end)];
            
            %plot interface
            
            figure(1)
            plot([yyy'; -flip(yyy')],-[xxx'; flip(xxx')]+xcm,'-')
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            cut = 3;
            zoomCut = 0.03;
            axis([-cut cut -cut cut])
            if zoom==1
                axis([-zoomCut zoomCut -zoomCut-1 -1+zoomCut])
            end
            if plotMesh==1
                hold on
                plot([y; -flip(y)],-[x; flip(x)]+xcm,'o')
                hold off
            end
            title(['IDshape=' num2str(IDshape(l)) ' IDdelta=' num2str(IDdelta(k)) ' Round=' num2str(findRound)])
            drawnow
        end

        if plotCurv==1
            figure(2)
            K = A.data(:,3)+A.data(:,4);
            %k = A.data(:,4);
            %plot curvature
            plot(x-xcm,A.data(:,3),'o-')
            hold on
            plot(x-xcm,A.data(:,4),'o-')
            plot(x-xcm,K,'o-')
            grid on
            xlabel('x')
            ylabel('k')
            drawnow
            hold off
        end

        if plotNormalVectors==1
            
            %compute the spline coeff
            [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric (x', y');

            [nx,ny] = normalToSPlines(bx,cx,dx,by,cy,dy);
            
            figure(3)
            subplot(2,1,1)
            
            %plot normal vector
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
        
        if plotForceFree==1
            figure(4)
            
            %compute the spline coeff
            [~, bx, cx, dx, ~, by, cy, dy] = spline_symmetric (x', y');

            
            K = A.data(:,3)+A.data(:,4);
            [nx,ny] = normalToSPlines(bx,cx,dx,by,cy,dy);
            
            Fdrop = int_axis_spline_symmetric(x',y',K.*nx);
            
            %plot force acting on droplet
            hold on
            plot(i,Fdrop,'ok-')
            %plot(x,K.*nx,'o-')
            grid on
            hold off
            xlabel('ite')
            ylabel('F')
            drawnow
        end

    end

    cd(here)

    if YesBreak==1

       V = V(1:i-1);  Area = Area(1:i-1); elong = elong(1:i-1); deltaNorm = deltaNorm(1:i-1);

    end
    
    if k==1 && l>1

            dt = num2str(dt);
            %if findRound==1
            if  numel(num2str(dt))==4
                resNext = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(1) '_elem=' elem '_dt=' dt '0000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
            elseif numel(num2str(dt))==5
                resNext = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(1) '_elem=' elem '_dt=' dt '000_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
            elseif numel(num2str(dt))==6
                resNext = [path 'Ca=' Ca '_lambda=' lambda '_IDshape=' num2str(IDshape(l)) '_IDdelta=' num2str(1) '_Round=' num2str(1) '_elem=' elem '_dt=' dt '00_loop=' totLoop '_RK=' RK '_CPUs=' CPUs ''];
            end
            cd(resNext)
            A = importdata('ParametersOptions.txt');
            dt = str2double(dt);
            Tadd = A.data(end)*dt*checkpoint;
            %Tadd = 1;
            cd(here)
            %elseif findRound>1
                %Tadd==0
                
            Tstart = Tstart+Tadd;
    end
    
    if k<m
        col = [colors(l) '--'];
    else
        col = colors(l);
    end
    
    %plot area variation norm and volume error
    R0 = nthroot(V(1)/4/pi*3,3);
    A0 = 4*pi*R0^2;
    %myNorm = sqrt((Area-A0)/(Area(1)-A0));
    myNorm = (Area-A0);
    if plotAreaNorm==1
    figure(5)
    plot(Tstart:dt*checkpoint*step:dt*(numel(myNorm)-1)*checkpoint+Tstart,myNorm(1:step:end),col)
    grid on
    ylabel('||A||')
    xlabel('t')
    hold on
    drawnow
    end

    if plotVolume==1
    figure(6)
    %plot colume error
    errV = (V-V(1))/V(1);
    plot(Tstart:dt*checkpoint*step:dt*(numel(myNorm)-1)*checkpoint+Tstart,errV(1:step:end))
    grid on
    ylabel('err_V')
    xlabel('t')
    hold on
    drawnow
    end

    %plot rear drop elongation
    figure(7)
    %plot(Tstart:dt*checkpoint*step:dt*(numel(myNorm)-1)*checkpoint+Tstart,(elong(1:step:end)-R)/R,col)
    %plot(Tstart:step*dt*checkpoint:dt*(numel(myNorm)-1)*checkpoint+Tstart,(elong(1:step:end)-R0)/R0)
    plot(Tstart:dt*checkpoint*step:dt*(numel(myNorm)-1)*checkpoint+Tstart,elong(1:step:end),col)
    %plot(Tstart:dt*checkpoint*step:dt*(numel(myNorm)-1)*checkpoint+Tstart,deltaNorm(1:step:end))
    grid on
    ylabel('rear elongation')
    xlabel('t')
    hold on
    drawnow
    
    if l==BreakShape
        break;
    end
    
    

    end
    
    %Tstart = Tstart+Tadd;

end











