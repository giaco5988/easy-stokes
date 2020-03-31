%post processing spherical motor

clear variables
close all

%parameters
TendUP = 1;
dtUp = 0.001;
res = 0;
thetaUp = linspace(pi/64,pi/16,20);
%thetaUp = thetaUp(1:70);
coeffDistUP = 2;
remeshUP = 1;
repUP = 3*ones(numel(thetaUp),1);
coeffRepUP = 1e6*ones(numel(thetaUp),1);
distRepUP = 0.1*ones(numel(thetaUp),1);
nPerLenghtUP = 5*ones(numel(thetaUp),1);
tolUp = 1e-3;
ODEup = 2;
betaUp = 1.4;
Hup = 1.85;
Lup = 10;
BNup = 16;
nCycle = 15;

%options
plotBubbleVel = 0;
plotMotorVel = 0;
plotDistanceWallBubble = 0;
plotDOF = 0;
plotVolume = 0;
plotLegend = 0;
dimensional = 0*ones(10,1); Tscale = 1e-2;   Lscale = 1e-6; Vscale = 1e-4;

cellLegend = cell(numel(repUP),1);

%initialize
maxVelBubble = zeros(numel(thetaUp),1);
maxVelMotor = zeros(numel(thetaUp),1);
distanceTravelled = zeros(numel(thetaUp),1);
timeExit = zeros(numel(thetaUp),1);
avgVel = zeros(numel(thetaUp),1);
Oxygen = zeros(numel(thetaUp),1);

countOut = 0;
for kkk = 1:numel(repUP)
    
    thisOut  = 0;
    
    display([num2str(kkk) ' of ' num2str(numel(repUP))])

    if res==1
        source = '~/Documents/MATLAB/droplet_simulations/results/';
    elseif res==0
        source = '~/Documents/MATLAB/droplet_simulations/server/';
    end

    BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorWithDiffusionSphericalBubble';
    here = pwd;

    if isempty(nCycle)
    
        filename = ['coneWithDiffusion_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUp) '_Hcc=' num2str(Hup) '_BN=' num2str(BNup) '_L=' num2str(Lup) '_nPerL=' num2str(nPerLenghtUP(kkk)) '_rep=' num2str(repUP(kkk)) '_coeffRep=' num2str(coeffRepUP(kkk)) '_distRep=' num2str(distRepUP(kkk)) '_theta=' num2str(thetaUp(kkk)) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];

        %upload
        load([source filename])

    else
   
        filename = ['coneWithDiffusion_nCycle=' num2str(nCycle) '_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUp) '_Hcc=' num2str(Hup) '_BN=' num2str(BNup) '_L=' num2str(Lup) '_nPerL=' num2str(nPerLenghtUP(kkk)) '_rep=' num2str(repUP(kkk)) '_coeffRep=' num2str(coeffRepUP(kkk)) '_distRep=' num2str(distRepUP(kkk)) '_theta=' num2str(thetaUp(kkk)) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];

        %upload
        load([source filename])

        %paste data
        for lll = 1:nCycle

            Tnow = manyT{lll};
            Ynow = manyY{lll};
            if lll==1
                Ynow = [Ynow(:,1:2) zeros(numel(Tnow),1) Ynow(:,3) zeros(numel(Tnow),1)];
                T = Tnow;
                Y = Ynow;
            else
                nBubble = (size(Ynow,2)-1)/2;

                if nBubble==2
                    T = [T; Tnow(2:end)];
                    Y = [Y; Ynow(2:end,:)];
                elseif nBubble==1
                    Ynow = [Ynow(2:end,1:2) zeros(numel(Tnow)-1,1) Ynow(2:end,3) zeros(numel(Tnow)-1,1)];
                    T = [T; Tnow(2:end)];
                    Y = [Y; Ynow];
                else
                    error('It is possible to have only one or two bubble')
                end
            end

        end
    
    end
    
    %upload
    load([source filename])
    
    cellLegend{kkk} = ['A=' num2str(coeffRepUP(kkk)) ' \delta=' num2str(distRepUP(kkk))];

    dist = zeros(size(Y,1),1);
    Vol = zeros(size(Y,1),1);
    distWallBubble = zeros(size(Y,1),1);
    countSnapshots = 1;
    tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp(kkk) linspace(0,1,100) linspace(-pi/2,pi/2,100)+thetaUp(kkk) linspace(0,pi,200)};
    for i = 1:size(Y,1)
        
        %find our if there is one or two bubble at this instant
        if Y(i,3)==0 && Y(i,5)==0
            bubbleNow = 1;
            tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp(kkk) linspace(0,1,100) linspace(-pi/2,pi/2,10)+thetaUp(kkk) linspace(0,pi,100)};
            PARAM.n = [100 10 100 10 100];
            PARAM.panels = [4 1];
        else
            bubbleNow = 2;
            tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUp(kkk) linspace(0,1,100) linspace(-pi/2,pi/2,10)+thetaUp(kkk) linspace(0,pi,100) linspace(0,pi,100)};
            PARAM.n = [100 10 100 10 100 100];
            PARAM.panels = [4 1 1];
        end

        %current location
        PARAM_now = PARAM;
        PARAM_now.xStart(1:4) = PARAM.xStart(1:4)+Y(i,1);
        PARAM_now.xEnd(1:4) = PARAM.xEnd(1:4)+Y(i,1);
        PARAM_now.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+Y(i,1);
        PARAM_now.x0_Circle(5) = Y(i,2);
        PARAM_now.rArc(5) = Y(i,4);
        if bubbleNow==2
            PARAM_now.x0_Circle(6) = Y(i,3);
            PARAM_now.rArc(6) = Y(i,5);
        end

        %compute curent shape
        [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);

        if plotDistanceWallBubble==1

            %distance wall-bubble
            [dist1,~,distWallBubble1] = panelDistance(x,y,5,3,PARAM_now);
            [dist2,~,distWallBubble2] = panelDistance(x,y,5,4,PARAM_now);
            dist(i) = min(dist1,dist2);
            distWallBubble(i) = min(distWallBubble1,distWallBubble2);

        end
        
        %compute bubble volume
        Vol(i) = 4/3*pi*PARAM_now.rArc(5)^3;

    end
    
    %compute oxygen consumption
    R1 = r;
    R2 = r+L*sin(thetaUp(kkk));
    h = L*cos(thetaUp(kkk));
    s = sqrt((R1-R2)^2+h^2);
    deltaT = T(end)-T(1);
    Oxygen(kkk) = pi*(R1+R2)*s*deltaT;
    
    %compute velocities
    D1 = finiteDifference1D(numel(T),[2 0],1);
    Umotor = D1*Y(:,1)./(D1*T);
    Ububble = D1*Y(:,2)./(D1*T);
    maxVelMotor(kkk) = max(abs(Umotor));
    maxVelBubble(kkk) = max(abs(Ububble));
    
    %compute distance travelled
    %deltaT = T(end)-T(1);
    distanceTravelled(kkk) = abs(Y(end,1))/Oxygen(kkk);

    if plotBubbleVel==1
        
        figure(1)
        
        if kkk>1
            hold on
        end

        plot(T,Ububble)
        xlabel('t')
        ylabel('U_b')
        title('Bubble velocity')
        grid on
        
    end

    if plotMotorVel==1
        
        figure(2)
        
        if kkk>1
            hold on
        end

        plot(T,Umotor)
        xlabel('t')
        ylabel('U_m')
        title('Motor velocity')
        grid on

    end

    if plotDistanceWallBubble==1
        
        figure(3)
        
        if kkk>1
            hold on
        end

        semilogy(T,dist,'-')
        xlabel('t')
        ylabel('d')
        %legend('min node distance','min shape dist','Location','Best')
        title('Distance wall-bubble')
        grid on

    end

    if plotDOF==1
        
        figure(4)
        
       if kkk>1
           hold on
       end

       plot(T,Y)
       grid on
       xlabel('t')
       title('Degrees of freedom')

    end
    
    %plot volume
    if plotVolume==1
        
       figure(5)
       
       if kkk>1
           hold on
       end
       
       semilogy(T,Vol)
       %plot(T,Vol)
       xlabel('t')
       ylabel('Volume')
       title('Bubble volume')
       grid on
        
    end
    
    timeExit(kkk) = T(end);
    avgVel(kkk) = 1/(T(end)-T(1))*trapz(T,Umotor);
    
    if thisOut==1
        
       avgThisOut(countOut) = avgVel(kkk);
       thetaOut(countOut) = thetaUp(kkk);
        
    end

    display(['Simulation time is ' num2str(simulationTime/60) ' min'])

end

if plotLegend==1

if plotBubbleVel==1
    figure(1)
    legend(cellLegend,'Location','Best')
end

if plotMotorVel==1
    figure(2)
    legend(cellLegend,'Location','Best')
end

if plotDistanceWallBubble==1
    figure(3)
    legend(cellLegend,'Location','Best')
end

if plotDOF==1
    figure(4)
    legend(cellLegend,'Location','Best')
end

if plotVolume==1
    figure(5)
    legend(cellLegend,'Location','Best')
end

end

figure
if dimensional==1
    plot(thetaUp,maxVelBubble*Vscale)
else
    plot(thetaUp,maxVelBubble)
end
grid on
xlabel('\theta')
if dimensional==1
    ylabel('max(Ub)[m/s]')
else
    ylabel('max(Ub)')
end
title('max bubble velocity')

figure
if dimensional==1
    plot(thetaUp,maxVelMotor*Vscale)
else
    plot(thetaUp,maxVelMotor)
end
grid on
xlabel('\theta')
if dimensional==1
    ylabel('max(Um)[m/s]')
else
    ylabel('max(Um)')
end
title('max motor velocity')

figure
if dimensional==1
    plot(thetaUp,timeExit*Tscale)
else
    plot(thetaUp,timeExit)
end
grid on
xlabel('\theta')
if dimensional==1
    ylabel('T_{exit}[s]')
else
    ylabel('T_{exit}')
end
title('Time to escape')

% figure
% if dimensional==1
%     plot(thetaUp,1./timeExit./Tscale)
% else
%     plot(thetaUp,1./timeExit)
% end
% grid on
% xlabel('\theta')
% if dimensional==1
%     ylabel('f_{ejection}[1/s]')
% else
%     ylabel('f_{ejection}')
% end
% title('Ejecion frequency')

figure
if dimensional==1
    plot(thetaUp,abs(avgVel)*Vscale)
else
    plot(thetaUp,abs(avgVel))
end
grid on
xlabel('\theta')
if dimensional==1
    ylabel('V_{avg}[m/s]')
else
    ylabel('V_{avg}')
end
title('Average motor velocity')
hold on
%plot(thetaOut,abs(avgThisOut),'or')

figure
if dimensional==1
    plot(thetaUp,distanceTravelled*Lscale)
else
    plot(thetaUp,distanceTravelled)
end
grid on
xlabel('\theta')
if dimensional==1
    ylabel('\Delta L[m]')
else
    ylabel('\Delta L')
end
title('Distance traveled over oxygen consumed')

% figure
% plot(thetaUp,abs(avgVel)./Oxygen)
% grid on
% xlabel('\theta')
% ylabel('V_{avg}/\Delta O_2')
% title('Average motor velocity over oxygen')

% figure
% plot(thetaUp,distanceTravelled./Oxygen)
% grid on
% xlabel('\theta')
% ylabel('\Delta L/\Delta O_2')
% title('Distance traveled per cycle over oxygen')






