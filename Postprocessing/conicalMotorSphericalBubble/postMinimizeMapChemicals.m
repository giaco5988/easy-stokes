%post processing spherical motor

clear variables
close all

set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%parameters
TendUP = 4;
dtUp = 1e-3;
res = 2;
%thetaUp = (2:2:20)/180*pi;
thetaUp = [10]/180*pi;
%thetaUp = 12/180*pi;
repUP = 2*ones(100,1);
coeffRepUP = 1e7*ones(100,1);
distRepUP = 0.1*ones(100,1);
nPerLenghtUP = 5*ones(100,1);
tolUp = 1e-3;
ODEup = 2;
betaUp = logspace(-1,1,64);
Hup = logspace(-3,1,64);
Lup = 10;
bubbleNucleatesUp = [];
nCycles = 0;
cycleType = 1;

%paramettric study
param2 = betaUp;   param2name = '\beta';
param1 = Hup;       param1name = 'H^{cc}';
thetaParam = thetaUp/pi*180;

%options
plotMapAVGvel = 1;  plotLinesOnMap = 1;
nanWhenIsIn = 1;

%parameters
Rstart = 1e-6;
gammaStart = 72*1e-3;
Astart = 1e-2;
Hstart = 0.1;
betaStart = 1;

%initialize
avgVel = zeros(numel(thetaUp),1);
Oxygen = zeros(numel(thetaUp),1);

%destination for saving plot
saveDest = '/Users/Giacomo/Documents/research_notes/APS_2017/movies/frames/';
namePlot = 'coneVelocityThree';
plotFig = 0;    savePlot = 0;   frameRate = 5;

color = get(gca,'ColorOrder');
close(figure(1))

maxR = zeros(numel(thetaUp),1);
maxGamma = zeros(numel(thetaUp),1);
maxV = zeros(numel(thetaUp),1);
for iii = 1:numel(thetaUp)

manyAVGvel = zeros(numel(param2),numel(param1));
for lll = 1:numel(param2)

countOut = 0;
thetaOut = [];
    
for kkk = 1:numel(param1)
    
    thisOut  = 0;
    thetaUpNow = thetaUp(iii);
    HccUPnow = Hup(kkk);
    betaUpNow = betaUp(lll);
    Lnow = Lup;
    
    display([num2str(kkk) ' of ' num2str(numel(repUP))])

    if res==1
        source = '~/Documents/MATLAB/droplet_simulations/results/';
    elseif res==0
        source = '~/Documents/MATLAB/droplet_simulations/server/';
    elseif res==2
        source = '~/Documents/MATLAB/droplet_simulations/results/microrocketWithDiffusionCambridge/';
    end

    BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorWithDiffusionSphericalBubble';
    here = pwd;

    %filename
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
    
    timeExit(kkk) = T(end);
    if numel(T)>2
        avgVel(kkk) = (Y(end,1)-Y(1,1))/(T(end)-T(1));
    else
        avgVel(kkk) = 0;
    end
    
    if thisOut==1 && nanWhenIsIn==1
        
        avgVel(kkk) = 0;
        
    end
    
    %store average vel
    manyAVGvel(lll,kkk) = avgVel(kkk);
    %totalDispCone(lll,kkk) = Y(end,1)-Y(1,1);
    totalDispCone(lll,kkk) = Y(end,1)-Y(1,1);

    display(['Simulation time is ' num2str(simulationTime/60) ' min'])

end

end

if plotMapAVGvel==1
    
    figure
    [~,h1] = contourf(log10(param1),log10(param2),abs(manyAVGvel),500);
    xlabel(['log' param1name])
    ylabel(['log' param2name])
    title('Average velocity')
    colorbar
    
    set(h1,'LineColor','none')
    
    %find global maximum for average velocity
    HmoveOnMap = @(var) (Hstart./(var(:,1)/Rstart).^2.*(var(:,2)/gammaStart))';
    betaMoveOnMap = @(var) (betaStart*(var(:,1)/Rstart)./(var(:,2)/gammaStart))';
    VaverageMAPdim = @(var) findClosestValueInterp(HmoveOnMap(var),betaMoveOnMap(var),Hup,betaUp,manyAVGvel).*(var(:,1)/Rstart)'./(var(:,2)/gammaStart)';
    maxRandGamma = fmincon(VaverageMAPdim,[Rstart gammaStart],[],[],[],[],[1e-1*Rstart 1e-1*gammaStart],[10*Rstart 10*gammaStart]);

    %compute new Hcc and beta and check if I am outside the map
    maxR(iii) = maxRandGamma(1);
    maxGamma(iii) = maxRandGamma(2);
    maxV(iii) = -VaverageMAPdim(maxRandGamma)*1e-4;
    HccNew = Hstart./(maxR(iii)/Rstart)^2*(maxGamma(iii)/gammaStart);
    betaNew = betaStart*(maxR(iii)/Rstart)/(maxGamma(iii)/gammaStart);
    display(['Maximum is obtained for Hcc=' num2str(HccNew) ' and beta=' num2str(betaNew)])
    if (HccNew>max(Hup) || HccNew<min(Hup)) || (betaNew>max(betaUp) || betaNew<min(betaUp))

        error('Maximum is ouside the map')

    end
    %hold on
    %plot(log10(HccNew),log10(betaNew),'.w','Markersize',40)
    
    %compute maximum varying only gamma or R
    nPlot = 500;
    RplotHere = Rstart*logspace(-1,1,nPlot);
    gammaPlotHere = gammaStart*logspace(-1,1,nPlot);
    %RplotHere = linspace(Rstart/10,Rstart*10,nPlot);
    %gammaPlotHere = gammaStart*linspace(gammaStart/10,gammaStart*10,nPlot);
    
    [temp1,ind1] = max(-VaverageMAPdim([RplotHere' gammaStart*ones(nPlot,1)])*1e-4);
    [temp2,ind2] = max(-VaverageMAPdim([Rstart*ones(nPlot,1) gammaPlotHere'])*1e-4);
    maxVelOnlyR(iii) = temp1;
    maxOnlyR(iii) = RplotHere(ind1);
    maxVelOnlyGamma(iii) = temp2;
    maxOnlyGamma(iii) = gammaPlotHere(ind2);
    
    if plotLinesOnMap==1
        
        nPlot = 500;
        
        %Rplot = linspace(1e-6*0.1,1e-5,nPlot);
        Rplot = Rstart*logspace(-1,1,nPlot);
        HplotLine = Hstart./(Rplot/Rstart).^2;
        betaPlotLine = betaStart*(Rplot/Rstart);
        [VlineR,fitR] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'k','Linewidth',3)
        
        %gammaPlot = linspace(72*1e-4,72*1e-2,nPlot);
        gammaPlot = gammaStart*logspace(-1,1,nPlot);
        HplotLine = Hstart*(gammaPlot/gammaStart);
        betaPlotLine = betaStart./(gammaPlot/gammaStart);
        [VlineGamma,fitGamma] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'r','Linewidth',3)
        
        %Aplot = linspace(Astart*0.01,Astart*100,100*nPlot);
        Aplot = Astart*logspace(-2,2,nPlot);
        HplotLine = Hstart./(Aplot/Astart);
        betaPlotLine = betaStart*ones(1,nPlot);
        [VlineFlux,fitFlux] = findClosestValueInterp(HplotLine,betaPlotLine,Hup,betaUp,abs(manyAVGvel));
        
        hold on
        plot(log10(HplotLine),log10(betaPlotLine),'w','Linewidth',3)
        
        hold on
        plot(log10(Hstart),log10(betaStart),'.k','Markersize',40)
        drawnow
        
        %plot dimensional velocity
        figure(2)
        if iii>1
            hold on
        end
        Rscale = Rplot./Rstart*1e-4;
        plot(Rplot,VlineR.*Rscale)
        grid on
        hold on
        plot(Rstart,fitR(Hstart,betaStart)*1e-4,'.k','Markersize',40)
        %plot(Rplot,Vline)
        xlabel('R[m]')
        ylabel('v[m/s]')
        
        %plot dimensional velocity
        figure(3)
        if iii>1
            hold on
        end
        gammaScale = gammaStart*1e-4./gammaPlot;
        plot(gammaPlot,VlineGamma.*gammaScale)
        grid on
        hold on
        plot(gammaStart,fitGamma(Hstart,betaStart)*1e-4,'.k','Markersize',40)
        %plot(Rplot,Vline)
        xlabel('\gamma [N/m]')
        ylabel('v[m/s]')
        
        %plot dimensional velocity
        figure(4)
        if iii>1
            hold on
        end
        fluxScale = Aplot/Astart*1e-4;
        loglog(Aplot,VlineFlux.*fluxScale)
        %plot(Aplot,VlineFlux)
        grid on
        hold on
        loglog(Astart,fitFlux(Hstart,betaStart)*1e-4,'.k','Markersize',40)
        %plot(Rplot,Vline)
        xlabel('A mol m^{-2} s^{-1}')
        ylabel('v[m/s]')
        
    end
    
end

end

if numel(thetaUp)>1
    
    figure
    plot(thetaParam,maxR)
    grid on
    xlabel('\theta')
    ylabel('R(max(V))[m]')
    
    figure
    plot(thetaParam,maxGamma)
    grid on
    xlabel('\theta')
    ylabel('\gamma (max(V))[N/m]')
    
    figure
    plot(thetaParam,maxV)
    grid on
    xlabel('\theta')
    ylabel('max(V)[m/s]')
    
    figure
    plot(thetaParam,maxOnlyR)
    grid on
    xlabel('\theta')
    ylabel('R(max(V)) \gamma=cosnt')
    
    figure
    plot(thetaParam,maxOnlyGamma)
    grid on
    xlabel('\theta')
    ylabel('\gamma(max(V)) R=const')

end

% for paper
figure
hold on
[ax,line1,line2] = plotyy(thetaParam,maxVelOnlyR,thetaParam,maxOnlyR);
grid on
xlabel('\theta')
ylabel(ax(1),'V_{max}')
ylabel(ax(2),'R_{max}')

%set(ax(1),'xlim')
line2.LineStyle = '-';
line2.LineStyle = '--';
set(line1,'Color','k')
set(line2,'Color','k')
set(line1,'Marker','o')
set(line2,'Marker','s')
set(line1,'MarkerSize',10)
set(line2,'MarkerSize',10)
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','k')
legend('Average velocity','Cone radius','Location','Best')
% set(ax(1),'xlim',[0 0.16])
% set(ax(2),'xlim',[0 0.16])
set(ax(1),'ylim',[4 5]*1e-4)
set(ax(2),'ylim',[1 2]*1e-6)
% %set(ax(1),'ylim',[0.049 0.06])
% %set(ax(2),'xlim',[0 26])
% %text(3.5,0.05,'\delta_{crit} prolate')
% set(ax(1),'XTick',0:0.04:0.16)
set(ax(1),'YTick',linspace(4,5,5)*1e-4)
set(ax(2),'YTick',linspace(1,2,5)*1e-6)























