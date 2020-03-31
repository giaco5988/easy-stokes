%post processing spherical motor

clear variables
close all

set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%parameters
TendUP = 1;
dtUp = 1e-4;
res = 0;
%thetaUp = linspace(pi/128,pi/16,20);
%thetaUp = linspace(pi/64,pi/6,10);
%thetaUp = [2 10 18];
%thetaUp = 2:2:20;
thetaUp = 2;
coeffDistUP = 2;
remeshUP = 1;
repUP = 2*ones(100,1);
coeffRepUP = 1e7*ones(100,1);
%coeffRepUP = [1e5 1e6 1e7 1e8];
distRepUP = 0.1*ones(100,1);
%distRepUP = 1e-1*ones(4,1);
nPerLenghtUP = 5*ones(100,1);
tolUp = 1e-3;
ODEup = 2;
%betaUp = logspace(-1,1,16);
%Hup = logspace(-3,1,16);
betaUp = [0.5 1 1.5];
%Hup = [1e-1 1 2];
%betaUp = 1;
Hup = 0.1;
Lup = 10;
bubbleNucleatesUp = [];
nCycles = 0;
cycleType = 1;

%paramettric study
%param1 = thetaUp;   param1name = '\theta';
%param1 = Lup;   param1name = 'L';
param1 = betaUp;   param1name = '\beta';
param2 = Hup;       param2name = 'H^{cc}';

%options
plotMapAVGvel = 0;  plotLinesOnMap = 1;
plotLegendParam2 = 0;
plotNoOut = 0;
nanWhenIsIn = 1;
plotBubbleVel = 0;  plotBubbleDisp = 1;
plotMotorVel = 0;   plotMotorDisp = 1;
plotLastBubblePos = 0;
plotDistanceWallBubble = 0;
plotDOF = 0;
plotVolume = 0; plotGrowthRate = 0;
plotRadius = 1; plotRadiusGrowth = 0;
plotLegend = 1;
dimensional = 02*ones(10,1); Tscale = 1e-2;   Lscale = 1e-6; Vscale = 1e-4;

cellLegend = cell(numel(param1),1);

%initialize
maxVelBubble = zeros(numel(thetaUp),1);
maxVelMotor = zeros(numel(thetaUp),1);
distanceTravelled = zeros(numel(thetaUp),1);
timeExit = zeros(numel(thetaUp),1);
avgVel = zeros(numel(thetaUp),1);
Oxygen = zeros(numel(thetaUp),1);

%destination for saving plot
saveDest = '/Users/Giacomo/Documents/research_notes/APS_2017/movies/frames/';
namePlot = 'coneVelocityThree';
plotFig = 0;    savePlot = 0;   frameRate = 5;

color = get(gca,'ColorOrder');

manyAVGvel = zeros(numel(param2),numel(param1));
for lll = 1:numel(param2)

countOut = 0;
thetaOut = [];

if plotLegendParam2==1
     cellLegend2{lll} = ['c_{crit}=' num2str(bubbleNucleatesUp(lll))];
end
    
for kkk = 1:numel(param1)
    
    thisOut  = 0;
    thetaUpNow = thetaUp/180*pi;
    HccUPnow = Hup(lll);
    betaUpNow = betaUp(kkk);
    Lnow = Lup;
    
    display([num2str(kkk) ' of ' num2str(numel(repUP))])

    if res==1
        source = '~/Documents/MATLAB/droplet_simulations/results/';
    elseif res==0
        source = '~/Documents/MATLAB/droplet_simulations/server/';
    end

    BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorWithDiffusionSphericalBubble';
    here = pwd;

    %filename
    %filename = ['coneWithDiffusion_nCycle=' num2str(nCycles) '_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUp) '_Hcc=' num2str(Hup) '_BN=' num2str(bubbleNucleatesUp(lll)) '_L=' num2str(Lup(kkk)) '_nPerL=' num2str(nPerLenghtUP(kkk)) '_rep=' num2str(repUP(kkk)) '_coeffRep=' num2str(coeffRepUP(kkk)) '_distRep=' num2str(distRepUP(kkk)) '_theta=' num2str(thetaUp(kkk)) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    filename = ['coneWithDiffusion_cycleType=' num2str(cycleType) '_nCycle=' num2str(nCycles) '_ODE=' num2str(ODEup) '_tol=' num2str(tolUp) '_beta=' num2str(betaUpNow) '_Hcc=' num2str(HccUPnow) '_BN=' num2str(bubbleNucleatesUp) '_L=' num2str(Lnow) '_nPerL=' num2str(nPerLenghtUP(kkk)) '_rep=' num2str(repUP(kkk)) '_coeffRep=' num2str(coeffRepUP(kkk)) '_distRep=' num2str(distRepUP(kkk)) '_theta=' num2str(thetaUpNow) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    
    %upload
    try
        load([source filename])
    catch
        Y = zeros(3,3);
        T = zeros(3,1);
        T(3) = 1;
        warning('File is empty')
    end
    
    %cellLegend{kkk} = ['A=' num2str(coeffRepUP(kkk)) ' \delta=' num2str(distRepUP(kkk))];
    if plotLegend==1
        cellLegend{kkk} = [param1name '=' num2str(param1(kkk))];
    end

    dist = zeros(size(Y,1),1);
    Vol = zeros(size(Y,1),1);
    distWallBubble = zeros(size(Y,1),1);
    countSnapshots = 1;
    tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUpNow linspace(0,1,100) linspace(-pi/2,pi/2,100)+thetaUpNow linspace(0,pi,200)};
    for i = 1:size(Y,1)
        
        %current location
        PARAM_now = PARAM;
        PARAM_now.xStart(1:4) = PARAM.xStart(1:4)+Y(i,1);
        PARAM_now.xEnd(1:4) = PARAM.xEnd(1:4)+Y(i,1);
        PARAM_now.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+Y(i,1);
        PARAM_now.x0_Circle(5) = Y(i,2);
        PARAM_now.rArc(5) = Y(i,3);

        %compute curent shape
        [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);

%         if PARAM_now.x0_Circle(5)>max(x{4})
%             display('Bubble is out')
%             countOut = countOut+1;
%             thetaOut(countOut) = thetaUp(kkk);
%             Lout(countOut) = Lup(kkk);
%             thisOut = 1;
%         end

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
    
    if T(end)==Tend && thisOut==0
       
        display('Simulation finishes before cycle ends')
        countOut = countOut+1;
        %param1(countOut) = thetaUp(kkk);
        %Lout(countOut) = Lup(kkk);
        thisOut = 1;
        
    end
    
    %compute oxygen consumption
    R1 = r;
    R2 = r+L*sin(thetaUpNow);
    h = L*cos(thetaUpNow);
    s = sqrt((R1-R2)^2+h^2);
    deltaT = T(end)-T(1);
    Oxygen(kkk) = pi*(R1+R2)*s*deltaT;
    
    %compute velocities
    if numel(T)>2
        D1 = finiteDifference1D(numel(T),[2 0],1);
    else
        D1 = [-1 1; -1 1];
    end
    Umotor = D1*Y(:,1)./(D1*T);
    Ububble = D1*Y(:,2)./(D1*T);
    if numel(Y)>3
        maxVelMotor(kkk) = max(abs(Umotor));
        maxVelBubble(kkk) = max(abs(Ububble));
    end
    
    %compute distance travelled
    distanceTravelled(kkk) = abs(Y(end,1));
    
    %last bubble position
    lastBubblePos(kkk) = Y(end,2)-Y(end,1);

    if plotBubbleVel==1
        
        figure(1)
        
        if kkk>1
            hold on
        end
        
%         if kkk==plotFig
%             
%             for indHere = 1:frameRate:numel(T)
%                 
%                 display([num2str(indHere) ' of ' num2str(numel(T))])
%         
%                 %subplot(2,1,2)
%                 hold on
%                 plot(T(1:indHere),Ububble(1:indHere),'Color',color(kkk,:))
%                 %hold on
%                 %plot(T(indHere),Umotor(indHere),'.','MarkerSize',40,'Color',color(kkk,:))
%                 grid on
%                 xlabel('t')
%                 ylabel('U_{bubble}')
%                 hold off
%                 axis([0 1 -5 120])
%                 title('Bubble Velocity')
%                 drawnow
% 
%                 if savePlot==1
%                     print('-dpng','-loose','-r100',[saveDest namePlot sprintf('%03d',round((indHere-1)/frameRate)) '.png'])
%                 end
%             
%             end
%             
%         else

        plot(T,Ububble)
        xlabel('t')
        ylabel('U_b')
        title('Bubble velocity')
        grid on
        
        %end
        
    end

    if plotMotorVel==1
        
        figure(2)
        
        if kkk>1
            hold on
        end
        
        if kkk==plotFig
            
            for indHere = 1:frameRate:numel(T)
                
                display([num2str(indHere) ' of ' num2str(numel(T))])
        
                %subplot(2,1,2)
                figure(2)
                hold on
                plot(T(1:indHere),Umotor(1:indHere),'Color',color(kkk,:))
                %hold on
                %plot(T(indHere),Umotor(indHere),'.','MarkerSize',40,'Color',color(kkk,:))
                grid on
                xlabel('t')
                ylabel('U_{cone}')
                hold off
                axis([0 1 -14 2])
                title('Cone Velocity')
                drawnow

                if savePlot==1
                    print('-dpng','-loose','-r100',[saveDest namePlot sprintf('%03d',round((indHere-1)/frameRate)) '.png'])
                end
            
            end
            
        else

            plot(T,Umotor)
            xlabel('t')
            ylabel('U_m')
            title('Cone velocity')
            grid on

        end

    end

    if plotDistanceWallBubble==1
        
        figure(3)
        
        if kkk>1
            hold on
        end

        %semilogy((0:i-2)*dt,distWallBubble)
        %hold on

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
    
    if plotVolume==1
        
        %compuet volume
        rBubble = Y(:,3);
        VolumeBubble = 4/3*pi*rBubble.^3;
        
        figure(5)
        
        if kkk>1
            hold on
        end

        plot(T,VolumeBubble)
        xlabel('t')
        ylabel('V')
        title('Bubble volume')
        grid on
        
    end
    
    if plotGrowthRate==1
        
        figure(6)
        
        %compuet volume
        rBubble = Y(:,3);
        VolumeBubble = 4/3*pi*rBubble.^3;
        growthRate = D1*VolumeBubble./(D1*T);
        
        if kkk>1
            hold on
        end

        plot(T,growthRate)
        xlabel('t')
        ylabel('dV/dt')
        title('Bubble growth rate')
        grid on
        
    end
    
    if plotRadius==1
        
        %compuet volume
        rBubble = Y(:,3);
        
        figure(7)
        
        if kkk>1
            hold on
        end
        
%         if kkk==plotFig
%             
%             for indHere = 1:frameRate:numel(T)
%                 
%                 display([num2str(indHere) ' of ' num2str(numel(T))])
%         
%                 %subplot(2,1,2)
%                 hold on
%                 plot(T(1:indHere),rBubble(1:indHere),'Color',color(kkk,:))
%                 grid on
%                 xlabel('t')
%                 ylabel('r_b')
%                 hold off
%                 axis([0 1 0 2])
%                 title('Bubble radius')
%                 drawnow
% 
%                 if savePlot==1
%                     print('-dpng','-loose','-r100',[saveDest namePlot sprintf('%03d',round((indHere-1)/frameRate)) '.png'])
%                 end
%             
%             end
%             
%         else

            plot(T,rBubble)
            xlabel('t')
            ylabel('r')
            title('Bubble radius')
            grid on
        
        %end
        
    end
    
    if plotRadiusGrowth==1
        
        figure(8)
        
        %compuet volume
        rBubble = Y(:,3);
        growthRateRadius = D1*rBubble./(D1*T);
        
        if kkk>1
            hold on
        end

        plot(T,growthRateRadius)
        xlabel('t')
        ylabel('dr/dt')
        title('Radius velocity')
        grid on
        
    end
    
    if plotBubbleDisp==1
        
        figure(11)
        
        if kkk>1
            hold on
        end

        plot(T,Y(:,2))
        %plot(Y(:,3),Y(:,2),'color',color(kkk,:))
        %xlabel('t')
        xlabel('r_b')
        ylabel('z_b')
        title('Bubble position')
        grid on
        
    end
    
    if plotMotorDisp==1
        
        figure(12)
        
        if kkk>1
            hold on
        end

        plot(T,Y(:,1))
        %plot(Y(:,3),Y(:,1),'color',color(kkk,:))
        %xlabel('t')
        xlabel('r_b')
        ylabel('z_c')
        title('Cone position')
        grid on
        
    end
    
    timeExit(kkk) = T(end);
    %avgVel(kkk) = 1/(T(i)-T(1))*trapz(T(1:i),Umotor(1:i));
    if numel(T)>2
        avgVel(kkk) = (Y(i,1)-Y(1,1))/(T(i)-T(1));
    else
        avgVel(kkk) = 0;
    end
    
    if thisOut==1 && nanWhenIsIn==1
        
        avgVel(kkk) = 0;
        
    end
    
    %store average vel
    manyAVGvel(lll,kkk) = avgVel(kkk);
    %totalDispCone(lll,kkk) = Y(end,1)-Y(1,1);
    totalDispCone(lll,kkk) = Y(i,1)-Y(1,1);
    timeToExit(lll,kkk) = T(i)-T(1);
    
%     if thisOut==1
%         
%        avgThisOut(countOut) = avgVel(kkk);
%        lastBubblePosOut(countOut) = lastBubblePos(kkk);
%        %break
%         
%     end

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

if plotGrowthRate==1
    figure(6)
    legend(cellLegend,'Location','Best')
end

if plotRadius==1
    figure(7)
    legend(cellLegend,'Location','Best')
end

if plotRadiusGrowth==1
    figure(8)
    legend(cellLegend,'Location','Best')
end

if plotBubbleDisp==1
    figure(11)
    legend(cellLegend,'Location','Best')
end

if plotMotorDisp==1
    figure(12)
    legend(cellLegend,'Location','Best')
end

end

    if plotNoOut==1
        [~,indMin] = min(abs(thetaUpNow-thetaOut(1)));
        indMin = indMin-1;
        thetaHere = thetaUp(1:indMin);
        avgVel = avgVel(1:indMin);
        lastBubblePos = lastBubblePos(1:indMin);
    else
        thetaHere = thetaUp;
    end

    figure(9)
    if lll>1
        hold on
    end
    if dimensional==1
        plot(param1,abs(avgVel)*Vscale,'x-')
    else
        plot(param1,abs(avgVel),'x-')
    end
    grid on
    xlabel(param1name)
    if dimensional==1
        ylabel('V_{avg}[m/s]')
    else
        ylabel('V_{avg}')
    end
%     if numel(bubbleNucleatesUp)==1
%         title(['Average motor velocity c_{crit}=' num2str(bubbleNucleates)])
%     else
        title('Average cone velocity')
    %end
    hold on
    if isempty(thetaOut)==0 && plotNoOut==0
        plot(param1,abs(avgThisOut),'or')
    end
    
    figure(10)
    if lll>1
        hold on
    end
    if dimensional==1
        plot(param1,abs(avgVel)*Vscale./Oxygen,'x-')
    else
        plot(param1,abs(avgVel)./Oxygen,'x-')
    end
    grid on
    xlabel(param1name)
    if dimensional==1
        ylabel('V_{avg}[m/s]')
    else
        ylabel('V_{avg}/O2')
    end
%     if numel(bubbleNucleatesUp)==1
%         title(['Average motor velocity c_{crit}=' num2str(bubbleNucleates)])
%     else
        title('Average per oxygen production')
    %end
    hold on
    if isempty(thetaOut)==0 && plotNoOut==0
        plot(param1,abs(avgThisOut),'or')
    end

%     figure(11)
%     if lll>1
%         hold on
%     end
%     plot(thetaHere,lastBubblePos./(L*cos(thetaHere)),'x-')
%     xlabel('\theta')
%     ylabel('$$ z_b/(L \cos \theta) $$','Interpreter','latex')
%     if numel(bubbleNucleatesUp)==1
%         title(['Last bubble position c_{crit}=' num2str(bubbleNucleates)])
%     else
%         title('Last bubble position')
%     end
%     grid on
%     hold on
%     if isempty(thetaOut)==0 && plotNoOut==0
%         plot(thetaOut,lastBubblePosOut./(L*cos(thetaOut)),'or')
%     end

    figure(13)
    if lll>1
        hold on
    end
    if dimensional==1
        plot(param1,1,'x-')
    else
        %plot(param1,abs(totalDispCone),'x-')
        [ax,line1,line2] = plotyy(param1,abs(totalDispCone),param1,1./timeToExit);
    end
    grid on
    xlabel(param1name)
    if dimensional==1
        ylabel('V_{avg}[m/s]')
    else
        %ylabel('|z_c|')
        ylabel(ax(1),'\Delta z_c')
        ylabel(ax(2),'1/T')
    end
    title('Final Cone Position')
    hold on
    if isempty(thetaOut)==0 && plotNoOut==0
        plot(param1,abs(avgThisOut),'or')
    end
    %axis([Hup(1) Hup(end) min(totalDispCone)+0.2*min(totalDispCone) 0])

end

if plotLegendParam2==1
    
    figure(7)
    legend(cellLegend2,'Location','Best')
    
    figure(8)
    legend(cellLegend2,'Location','Best')
    
end

if plotMapAVGvel==1
    
    figure
    [~,h1] = contourf(log10(param1),log10(param2),abs(manyAVGvel),500);
    xlabel(['log' param1name])
    ylabel(['log' param2name])
    title('Average velocity')
    colorbar
    
    set(h1,'LineColor','none')
    
    if plotLinesOnMap==1
        
        nPlot = 500;
        Hstart = 0.1;
        betaStart = 1;
        
        Rstart = 1e-6;
        %Rplot = linspace(1e-6*0.1,1e-5,nPlot);
        Rplot = Rstart*logspace(-1,1,nPlot);
        HplotLine = Hstart./(Rplot/Rstart).^2;
        betaPlotLine = betaStart*(Rplot/Rstart);
        [VlineR,fitR] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'k','Linewidth',3)
        
        gammaStart = 72*1e-3;
        %gammaPlot = linspace(72*1e-4,72*1e-2,nPlot);
        gammaPlot = gammaStart*logspace(-1,1,nPlot);
        HplotLine = Hstart*(gammaPlot/gammaStart);
        betaPlotLine = betaStart./(gammaPlot/gammaStart);
        [VlineGamma,fitGamma] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'r','Linewidth',3)
        
        Astart = 1e-2;
        %Aplot = linspace(Astart*0.01,Astart*100,100*nPlot);
        Aplot = Astart*logspace(-2,2,nPlot);
        HplotLine = Hstart./(Aplot/Astart);
        betaPlotLine = betaStart*ones(1,nPlot);
        [VlineFlux,fitFlux] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'w','Linewidth',3)
        
        hold on
        plot(log10(Hstart),log10(betaStart),'.k','Markersize',40)
        
        %plot dimensional velocity
%         figure(30)
%         hold on
%         Rscale = Rplot./Rstart*1e-4;
%         plot(Rplot,VlineR.*Rscale)
%         grid on
%         hold on
%         plot(Rstart,fitR(Hstart,betaStart)*1e-4,'.k','Markersize',40)
%         %plot(Rplot,Vline)
%         xlabel('R[m]')
%         ylabel('v[m/s]')
%         
%         %plot dimensional velocity
%         figure(31)
%         hold on
%         gammaScale = gammaStart*1e-4./gammaPlot;
%         plot(gammaPlot,VlineGamma.*gammaScale)
%         grid on
%         hold on
%         plot(gammaStart,fitGamma(Hstart,betaStart)*1e-4,'.k','Markersize',40)
%         %plot(Rplot,Vline)
%         xlabel('\gamma [N/m]')
%         ylabel('v[m/s]')
%         
%         %plot dimensional velocity
%         figure(32)
%         hold on
%         fluxScale = Aplot/Astart*1e-4;
%         loglog(Aplot,VlineFlux.*fluxScale)
%         %plot(Aplot,VlineFlux)
%         grid on
%         hold on
%         loglog(Astart,fitFlux(Hstart,betaStart)*1e-4,'.k','Markersize',40)
%         %plot(Rplot,Vline)
%         xlabel('A mol m^{-2} s^{-1}')
%         ylabel('v[m/s]')
        
        %plot dimensional velocity
        figure(30)
        hold on
        Rscale = Rplot./Rstart*1e-4;
        [axR,line1R,line2R] = plotyy(Rplot,VlineR,Rplot,VlineR.*Rscale);
        %plot(Rplot,VlineR.*Rscale)
        grid on
        hold on
        %[axRsim,line1Rsim,line2Rsim] = plotyy(Rstart,fitR(Hstart,betaStart),Rstart,fitR(Hstart,betaStart)*1e-4);
%         set(line1Rsim,'Marker','.')
%         set(line1Rsim,'Color','k')
%         set(line1Rsim,'MarkerSize',40)
%         set(line2Rsim,'Marker','.')
%         set(line2Rsim,'Color','k')
%         set(line2Rsim,'MarkerSize',40)
        
        %plot(Rplot,Vline)
        xlabel('R[m]')
        %ylabel('v[m/s]')
        ylabel(axR(1),'U')
        ylabel(axR(2),'U_{dim}')
        set(axR(1),'YTick',0:2:10)
        set(axR(2),'YTick',0:1e-4:5e-4)
        
        line2R.LineStyle = '--';
        set(line1R,'Color','k')
        %set(line1,'Marker','o')
        set(line2R,'Color','k')
        set(axR(1),'ycolor','k')
        set(axR(2),'ycolor','k')
        
        %plot dimensional velocity
        figure(31)
        hold on
        gammaScale = gammaStart*1e-4./gammaPlot;
        %plot(gammaPlot,VlineGamma.*gammaScale)
        [axGamma,line1Gamma,line2Gamma] = plotyy(gammaPlot,VlineGamma,gammaPlot,VlineGamma.*gammaScale);
        grid on
        hold on
        %plotyy(gammaStart,fitGamma(Hstart,betaStart)*1e-4,gammaStart,fitGamma(Hstart,betaStart))
        %plot(gammaStart,fitGamma(Hstart,betaStart),'.k','Markersize',40)
        %plot(Rplot,Vline)
        xlabel('\gamma [N/m]')
        %ylabel('v[m/s]')
        ylabel(axGamma(1),'U')
        ylabel(axGamma(2),'U_{dim}')
        
        line2Gamma.LineStyle = '--';
        set(line1Gamma,'Color','k')
        %set(line1,'Marker','o')
        set(line2Gamma,'Color','k')
        set(axGamma(1),'ycolor','k')
        set(axGamma(2),'ycolor','k')
        set(axGamma(1),'YTick',0:2:10)
        set(axGamma(2),'YTick',0:1e-4:5e-4)
        
    end
    
%     figure
%     surf(param1,param2,abs(manyAVGvel));
%     xlabel(param1name)
%     ylabel(param2name)
%     zlabel('|V_{avg}|')
%     title('Average velocity')
%     %colorbar
    
%     figure
%     [~,h1] = contourf(param1,param2,abs(manyAVGvel)./repmat(Oxygen,1,numel(param2))',500);
%     xlabel(param1name)
%     ylabel(param2name)
%     title('Average velocity/O2')
%     colorbar
%     
%     set(h1,'LineColor','none')
    
end

% figure(13)
% line2.LineStyle = '--';
% set(line1,'Color','k')
% %set(line1,'Marker','o')
% set(line2,'Color','k')
% set(ax(1),'ycolor','k')
% set(ax(2),'ycolor','k')
% %set(ax(1),'xlim',[0 0.16])
% %set(ax(2),'xlim',[0 0.16])
% set(ax(1),'ylim',[0.6 1.8])
% set(ax(2),'ylim',[2 8])
% %set(ax(1),'ylim',[0.049 0.06])
% %set(ax(2),'xlim',[0 26])
% %text(3.5,0.05,'\delta_{crit} prolate')
% %set(ax(1),'XTick',0:0.04:0.16)
% set(ax(1),'YTick',0.6:0.2:1.8)
% set(ax(2),'YTick',2:8)























