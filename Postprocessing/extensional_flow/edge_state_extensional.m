%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variables
Ca = '0.040000';
lambda = '1.000000';

% options
plotShape = 0;  zoom = 0;   plotMesh = 0;
plotCurv = 0; 
plotNormalVectors = 0;
plotVolume = 0;
plotAreaNorm = 0;
step = 1;
plotRES = 0;

% options serie computation
ComputeLegendre = 0;
symmetric = 1;
modes = 3;

colors= 'brgkcm';
%colors= 'bbbbbb';
%colors= 'rrrrrr';
colors = repmat(colors,1,1000);

%path
here = pwd;
%path = '~/Documents/C++/edge_state_extensional/results/';
%path = '~/Documents/C++/edge_state_extensional/resultsServer2/';
path = '~/Documents/C++/edge_state_extensional/resultsCaSmall/';

%ID
IDshape = 2:5;    BreakShape = Inf;
IDdelta = 1:16;
Round = 1;
CheckRound = 1;
numIDdelta = 16;
only2 = 0;
elem = '200';
dt = 0.01;
totLoop = '25000';
RK = '2';
CPUs = '1';

%get steady state solution form newton
load('~/Documents/MATLAB/Postprocessing/extensional_flow/CaExt.mat')
%load('~/Documents/MATLAB/Postprocessing/extensional_flow/DExt.mat')
load('~/Documents/MATLAB/Postprocessing/extensional_flow/manyA.mat')
%load('~/Documents/MATLAB/Postprocessing/extensional_flow/manyB.mat')

[~,indSide] = min(abs(str2double(Ca)-manyCa));
xSteadyStable = max(manyA(:,indSide));

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
    Dellipse = zeros(loop,1);
    V = zeros(loop,1);
    elong = zeros(loop,1);
    deltaNorm = zeros(loop,1);
    fMode = zeros(loop,modes);
    errSerie = zeros(loop,1);
    
    %define time scale
    Tscale = 1/str2double(Ca);
    WhenRes = Tscale/dt/checkpoint;

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
        
        %place drop in the center
        %x = x-xcm;

        %compute and save area
        Area(i) = surf_gauss_vect(x',y');
        
        %ellipsicity
        major = max(x-xcm); minor = max(y);
        Dellipse(i) = (major-minor)/(major+minor);

        %compute and save volume
        V(i) = axis_int_gauss_vect(x',y');
        
        %compute legendre serie
        if ComputeLegendre==1;
            [fMode(i,:),errSerie(i)] = LegendreSerie(x-xcm,y,modes,symmetric,V(i));
        end

        %compute lateral elongation compared to stable shape
        elong(i) = x(1) - xcm - xSteadyStable;
        
        %compute norm of perturbation
        R0 = nthroot(V(i)/4/pi*3,3);
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
            %subplot(3,1,1)
            plot([xxx'; flip(xxx')],[yyy'; -flip(yyy')],'-')
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            cut = 6;
            zoomCut = 0.03;
            axis([-cut cut -cut cut])
            if zoom==1
                axis([-zoomCut zoomCut -zoomCut-1 -1+zoomCut])
            end
            if plotMesh==1
                hold on
                plot([x; flip(x)],[y; -flip(y)],'o')
                hold off
            end
            title(['IDshape=' num2str(IDshape(l)) ' IDdelta=' num2str(IDdelta(k)) ' Round=' num2str(findRound)])
            drawnow
            
%             theta = atan(yyy./xxx);
%             theta = theta + pi*(theta<0);
%             
%             figure(20)
%             %subplot(3,1,2)
%             plot(theta,xxx,'k')
%             xlabel('\theta')
%             ylabel('x')
%             grid on
%             drawnow
%             
%             figure(21)
%             %subplot(3,1,3)
%             plot(theta,yyy,'k')
%             xlabel('\theta')
%             ylabel('y')
%             grid on
%             drawnow
%             
%             figure(22)
%             %subplot(3,1,3)
%             plot(theta,sqrt(xxx.^2+yyy.^2),'k')
%             xlabel('\theta')
%             ylabel('y')
%             grid on
%             drawnow
            
        end

        if plotCurv==1
            figure(2)
            K = A.data(:,3)+A.data(:,4);
            %k = A.data(:,4);
            %plot curvature
            plot(x-xcm,K,'o-')
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

       V = V(1:i-1);  Area = Area(1:i-1); elong = elong(1:i-1); deltaNorm = deltaNorm(1:i-1);
       fMode = fMode(1:i-1,:);

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
        colPhase = [colors(l) '--'];
    else
        col = [colors(l) '-'];
        colPhase = [colors(l) '-'];
    end
    
    %plot area variation norm and volume error
    R0 = nthroot(V(1)/4/pi*3,3);
    A0 = 4*pi*R0^2;
    %myNorm = sqrt((Area-A0)/(Area(1)-A0));
    myNorm = (Area-A0);
    %myNorm = myNorm(1:i-1);
    if plotAreaNorm==1
        figure(5)
        xCoord = 0:numel(myNorm)-1;
        plot(Tstart+dt*checkpoint*step*xCoord,myNorm(1:step:end),col)
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
        xCoord = 0:numel(myNorm)-1;
        plot(Tstart+dt*checkpoint*step*xCoord,errV(1:step:end))
        grid on
        ylabel('err_V')
        xlabel('t')
        hold on
        drawnow
    end

    %plot side drop elongation
    figure(7)
    xCoord = 0:numel(myNorm)-1;
    plot(Tstart+dt*checkpoint*step*xCoord,elong(1:step:end),col)
    grid on
    ylabel('side elongation')
    title(['Edge tracking Ca=' Ca])
    xlabel('t')
    hold on
    drawnow
    
    %plot phase space
    if ComputeLegendre==1
        figure(9)
        %colHere = [col ''];
        if modes==3
            plot(fMode(:,2),fMode(:,3),colPhase)
            ylabel('f_2')
            xlabel('f_1')
        elseif modes==4
            plot3(fMode(:,2),fMode(:,3),fMode(:,4),colPhase)
            ylabel('f_2')
            xlabel('f_1')
            zlabel('f_3')
        end
        title(['Phase space Ca=' Ca])
        grid on
        hold on
        drawnow
    end
    
    %plot residuals
    if plotRES==1
        figure(8)
        Dellipse = Dellipse(1:i-1);
        res = sqrt((Dellipse(WhenRes+1:end)-Dellipse(1:end-WhenRes)).^2);
        semilogy((1:numel(Dellipse)-WhenRes)*checkpoint,res,col)
        xlabel('iteration')
        ylabel('Res_D')
        grid on
        hold on
        drawnow
    end
    
    if l==BreakShape
        break;
    end
    
    

    end
    
    %Tstart = Tstart+Tadd;

end











