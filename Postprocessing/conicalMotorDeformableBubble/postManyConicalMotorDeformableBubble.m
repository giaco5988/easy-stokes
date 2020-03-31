%post processing of a conical motor with a defromable bubble

clear variables
close all

%parameters
r = 0;
L = 10;
res = 0;
CaUP = 1e-2;
%CaUP = [1e-2 5e-2 1e-1];
CaUP = [1e-3 5e-3 1e-2 5e-2 1e-1];
thetaUp = 1;
%thetaUp = [0.25:0.25:1.75 2 3 4];
%thetaUp = [0.25 1 2];
%dtUp = [5e-2 2e-2 1e-2 1e-2];
dtUp = 2e-2;
nElemUP = 124;
TendUP = [10000 5000*ones(1,5)];
%TendUP = 5000;
coeffRepUP = -10;
repOn = 0.2;
alphaUP = 0.5;
%alphaUP = r+L/2*sin(thetaUp)-2*repOn;
ODE = 0;
repTypeUP = 8;
massFlux = 1;
%x0UP = 1:4;
x0UP = 5;

%parametric study on this one
%param1 = thetaUp;   param1name = '\theta';
param1 = CaUP;   param1name = '\rm{Ca}';
%param1 = x0UP;   param1name = 'z_0';

%and this one
%param2 = CaUP;   param2name = '\rm{Ca}';
param2 = thetaUp;   param2name = '\theta';

%options
plotMotorUp = 0;    plotIteUp = 1;
threeD = 0; theta3D = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.3];
plotLegend = 1;
plotConeVel = 1;
stopSphere = 1; stopValue = 2e-3+eps;
plotBubbleArea = 1;
plotBubbleXcm = 1;
plotDimRadius = 0;

%scales
gamma = 7.2*1e-2;       % with surfactants
Radius = 1e-6;          % microrocket radius
fluxChemical = 1e-2;    % chemical flux from Ref. Manjare
mu = 1e-3;              % water viscosity
Lscale = Radius;
Vscale = gamma/mu;
Tscale = Radius*mu/gamma;

maxAvgVel = zeros(numel(param2),1);
maxAvgVel1 = zeros(numel(param2),1);
maxAvgVel2 = zeros(numel(param2),1);
thetaMaxAvgVel = zeros(numel(param2),1);
thetaMaxAvgVel1 = zeros(numel(param2),1);
thetaMaxAvgVel2 = zeros(numel(param2),1);
for kkk = 1:numel(param2)

%initialize
Texit = zeros(numel(param1),1);
TbubblePartiallyOut = zeros(numel(param1),1);
ZCbubblePartiallyOut = zeros(numel(param1),1);
vAverage = zeros(numel(param1),1);
finalPos = zeros(numel(param1),1);
cellLegend = cell(numel(param1),1);
for lll = 1:numel(param1)
    
    thetaUPnow = thetaUp;
    alphaUPnow = alphaUP;
    CaUPnow = CaUP(lll);
    TendUPnow = TendUP(lll);
    nElemUPnow = nElemUP;
    x0UPnow = x0UP;
    dtUPnow = dtUp;
    cellLegend(lll) = {[param1name '=' num2str(param1(lll))]};

    %filename
    if ODE==0
        %filename = ['coneBubbleDeformable_x0=5_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUP) '_L=10_nStart=' num2str(nElemUP) '_rep=3_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUPnow) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
        filename = ['coneBubbleDeformable_massFlux=' num2str(massFlux) '_ODE=' num2str(ODE) '_x0=' num2str(x0UPnow) '_alpha=' num2str(alphaUPnow) '_Ca=' num2str(CaUPnow) '_L=10_nStart=' num2str(nElemUPnow) '_rep=' num2str(repTypeUP) '_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUPnow) '_Tend=' num2str(TendUPnow) '_dt=' num2str(dtUPnow) '.mat'];
    else
        filename = ['coneBubbleDeformable_ODE=' num2str(ODE) '_x0=5_alpha=' num2str(alphaUP) '_Ca=' num2str(CaUPnow) '_L=10_nStart=' num2str(nElemUP) '_rep=3_coeffRep=' num2str(coeffRepUP) '_distRep=0.2_theta=' num2str(thetaUPnow) '_Tend=' num2str(TendUP) '_dt=' num2str(dtUp) '.mat'];
    end

    %source
    if res==1
       source = '~/Documents/MATLAB/droplet_simulations/results/';
    elseif res==0
       source = '~/Documents/MATLAB/droplet_simulations/server/';
    elseif res==2
       source = '~/Documents/MATLAB/test/coneWithDeformableBubble/resTest/';
    end

    BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/conicalMotorDeformablelBubble';
    here = pwd;

    %load file
    load([source filename])
    PARAM.ellipseShape = [0 0 0 0 1];

    %number
    if T(end)~=0
        loopUp = numel(T);
    else
        loopUp = find(T==0,2);
        loopUp = loopUp(2)-1;
    end

    %loop
    bubbleIsOut = 0;
    bubbleIsConfined = 0;
    bubbleIsExiting = 0;
    Ucone = zeros(numel(T),1);
    manyPosMotor = zeros(numel(T),1);
    Area = zeros(numel(T),1);
    Volume = zeros(numel(T),1);
    xcm = zeros(numel(T),1);
    tParametricPost = {linspace(0,1,100) linspace(pi/2,3*pi/2,10)+thetaUPnow/180*pi linspace(0,1,100) linspace(-pi/2,pi/2,100)+thetaUPnow/180*pi linspace(0,pi,200)};
    for i = 1:loopUp

       display([num2str(i) ' of ' num2str(numel(T))])

       %current shape and motor position
       Yhere = Y{i};
       xBubble = Yhere(1:2:end-2);
       yBubble = Yhere(2:2:end-1);
       posMotor = Yhere(end);

       %current location
       PARAM_now = PARAM;
       PARAM_now.xStart(1:4) = PARAM.xStart(1:4)+posMotor;
       PARAM_now.xEnd(1:4) = PARAM.xEnd(1:4)+posMotor;
       PARAM_now.x0_Circle(1:4) = PARAM.x0_Circle(1:4)+posMotor;

       %compute curent shape
       [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);
       x{5} = xBubble';
       y{5} = yBubble';

       %motor velocity
       if ODE==0
        Vhere = V{i};
        Ucone(i) = Vhere(end);
       end
       manyPosMotor(i) = Yhere(end);

       %bubble surface area
       Area(i) = surf_gauss_vect(xBubble',yBubble');
       Volume(i) = axis_int_gauss_vect(xBubble',yBubble');

       %bubble center of mass
       xcm(i) = center_mass(xBubble',yBubble');
       %xcm(i) = centerOfMassBlockAxis(x,y,2,PARAM);
       
       if bubbleIsOut==0 && max(xBubble)>max(x{4})
           
          indBubbleIsPartiallyOut = i;
          bubbleIsOut = 1;
           
       end
       
       if bubbleIsConfined==0 && max(xBubble)<max(x{4})
           
           aBubble = max(xBubble)-min(xBubble);
           bBubble = 2*max(yBubble);
           Dbubble = (aBubble-bBubble)/(aBubble+bBubble);

           if Dbubble>stopValue*10
           
            indBubbleIsConfined = i;
            bubbleIsConfined = 1;
          
           end
           
       end
       
       if bubbleIsExiting==0 && max(xBubble)>max(x{4})
           
            indBubbleIsExiting = i;
            bubbleIsExiting = 1;
           
       end

       %plot shape
       if plotMotorUp==1 && sum(i==(1:plotIteUp:numel(T)))

            %block coordinates
            figure(1)
            PARAM_now.eliminateBlockPostprocessing = [0 0];
            plotGeometryStokesV3(x,y,threeD,theta3D,color3D,transp,[1 1],PARAM_now)
            axis([-5 15 -7 7])
            grid on
            hold off
            axis equal
            axis off
            title(['t=' num2str(T(i))])
            drawnow

       end

       if stopSphere==1 && max(xBubble)>max(x{4})

           aBubble = max(xBubble)-min(xBubble);
           bBubble = 2*max(yBubble);
           Dbubble = (aBubble-bBubble)/(aBubble+bBubble);

           if Dbubble<stopValue

               disp('Bubble is out and almost spherical')
               break

           end

       end

    end

    %derivative
    D1 = finiteDifference1D(i,[2 0],1);
    
    if ODE~=0
        Ucone = (D1*manyPosMotor(1:i))./(D1*T(1:i));
    end
    
    %compute time to exit and final cone position
    %TbubblePartiallyOut(lll) = T(indBubbleIsPartiallyOut)-T(indBubbleIsConfined);
    %ZCbubblePartiallyOut(lll) = manyPosMotor(indBubbleIsPartiallyOut)-manyPosMotor(indBubbleIsConfined);
    
    %compute time to exit and final cone position
    Texit(lll) = T(i);
    finalPos(lll) = manyPosMotor(i)-manyPosMotor(1);
    
    %compute average velocity
    %vAverage(lll) = 1/(T(i)-T(1))*trapz(T(1:i),Ucone(1:i));
    vAverage(lll) = (manyPosMotor(i)-manyPosMotor(1))/(T(i)-T(1));

    %plot bubble surface area
    if plotBubbleArea==1

        figure(2)
        if lll>1
            hold on
        end
        %subplot(2,1,1)
        plot(T(1:i),Area(1:i))
        grid on
        xyLabelTex('t','A')
        title('Surface area')
        drawnow
        
        %excess surface area
        figure(9)
        if lll>1
            hold on
        end
        %subplot(2,1,1)
        Rsphere = nthroot(3/4/pi*Volume(1:i),3);
        excessArea = Area(1:i)-4*pi*Rsphere.^2;
        plot(T(1:i),excessArea)
        grid on
        xyLabelTex('t','\Delta A')
        title('Excess surface area')
        drawnow
        
        figure(31)
        if lll>1
            hold on
        end
        [~,indMax] = max(abs(excessArea(1:i)));
        plot(T(1:i)/T(indMax),excessArea(1:i)/excessArea(indMax))
        grid on
        xyLabelTex('t/t^*_A','\Delta A/\Delta A^*_A')
        %xlabel('t')
        %ylabel('U_M')
        title('Rescaled excess surface area')
        drawnow

        dAdt = (D1*Area(1:i))./(D1*T(1:i));

        %subplot(2,1,2)
        figure(3)
        if lll>1
            hold on
        end
        plot(T(1:i),dAdt(1:i))
        grid on
        xyLabelTex('t','dA/dt')
        %xlabel('t')
        %ylabel('dA/dt')
        title('Surface area variation')
        drawnow

    end
    
    %identify recoling phase
    [maxExcessArea,startToRecoil] = max(excessArea);
    figure(9)
    hold on
    plot(T(startToRecoil),maxExcessArea,'ok')
    plotMaxExcessArea(lll) = maxExcessArea;
    
    %compute time to exit and final cone position
    TbubblePartiallyOut(lll) = T(startToRecoil)-T(1);
    ZCbubblePartiallyOut(lll) = manyPosMotor(startToRecoil)-manyPosMotor(1);

    %plot bubble center of mass
    if plotBubbleXcm==1

        figure(4)
        if lll>1
            hold on
        end
        %subplot(2,1,1)
        plot(T(1:i),xcm(1:i))
        grid on
        xyLabelTex('t','z_b^{\rm{cm}}')
        %xlabel('t')
        %ylabel('z_b^{\text{cm}}')
        title('Bubble position')
        drawnow
        
        figure(32)
        if lll>1
            hold on
        end
        %subplot(2,1,1)
        plot(T(1:i)/T(i),xcm(1:i)/xcm(i))
        grid on
        xyLabelTex('t/t_{final}','z_b/z_b^{final}')
        %xlabel('t')
        %ylabel('z_b^{\text{cm}}')
        title('Bubble position rescaled')
        drawnow

        dXdt = (D1*xcm(1:i))./(D1*T(1:i));

        %subplot(2,1,2)
        figure(5)
        if lll>1
            hold on
        end
        plot(T(1:i),dXdt(1:i))
        grid on
        xyLabelTex('t','\dot{z}_b')
        %xlabel('t')
        %ylabel('U_b')
        title('Bubble velocity')
        drawnow
        
        figure(30)
        if lll>1
            hold on
        end
        [~,indMax] = max(abs(dXdt(1:i)));
        plot(T(1:i)/T(indMax),dXdt(1:i)/dXdt(indMax))
        grid on
        xyLabelTex('t/t^*_B','\dot{z}_b/\dot{z}_b^*')
        %xlabel('t')
        %ylabel('U_M')
        title('Rescaled bubble velocity')
        drawnow

    end

    %plot motor velocity
    if plotConeVel==1

        figure(6)
        if lll>1
            hold on
        end
        plot(T(1:i),Ucone(1:i))
        grid on
        xyLabelTex('t','\dot{z}_c')
        %xlabel('t')
        %ylabel('U_M')
        title('Cone velocity')
        drawnow
        
        % fit power law to cone velocity in the migration phase
        deltaTmigration(lll) = T(indBubbleIsExiting)-T(indBubbleIsConfined);
        figure(40)
        loglog(T(indBubbleIsConfined:indBubbleIsExiting)-T(indBubbleIsConfined),abs(Ucone(indBubbleIsConfined:indBubbleIsExiting)-Ucone(indBubbleIsConfined)))
        hold on
        grid on
        fitPowerLaw = fit(T(indBubbleIsConfined+1:indBubbleIsExiting)-T(indBubbleIsConfined),abs(Ucone(indBubbleIsConfined+1:indBubbleIsExiting)-Ucone(indBubbleIsConfined)),'power1');
        coeffPL(lll) = fitPowerLaw.a;
        exponentPL(lll) = fitPowerLaw.b;
        
        figure(29)
        if lll>1
            hold on
        end
        [~,indMax] = max(abs(Ucone(1:i)));
        plot(T(1:i)/T(indMax),Ucone(1:i)/Ucone(indMax))
        grid on
        xyLabelTex('t/t^*_C','\dot{z}_c/\dot{z}_c^*')
        %xlabel('t')
        %ylabel('U_M')
        title('Rescaled cone velocity')
        drawnow
        
        figure(7)
        if lll>1
            hold on
        end
        plot(T(1:i),manyPosMotor(1:i))
        grid on
        xyLabelTex('t','z_c')
        %xlabel('t')
        %ylabel('U_M')
        title('Cone Position')
        drawnow
        
        figure(8)
        if lll>1
            hold on
        end
        plot(T(1:i)/T(i),manyPosMotor(1:i)/manyPosMotor(i))
        grid on
        xyLabelTex('t/t_{final}','z_c/z_c^{final}')
        %xlabel('t')
        %ylabel('U_M')
        title('Cone Position Rescaled')
        drawnow
        
        deltaT1(lll) = T(indBubbleIsConfined)-T(1);
        deltaT2(lll) = T(startToRecoil)-T(indBubbleIsConfined);
        deltaT3(lll) = T(i)-T(startToRecoil);
        deltaZ1(lll) = abs(manyPosMotor(indBubbleIsConfined)-manyPosMotor(1));
        deltaZ2(lll) = abs(manyPosMotor(startToRecoil)-manyPosMotor(indBubbleIsConfined));
        deltaZ3(lll) = abs(manyPosMotor(i)-manyPosMotor(startToRecoil));
        Vavg_1(lll) = deltaZ1(lll)./deltaT1(lll);
        Vavg_2(lll) = deltaZ2(lll)./deltaT2(lll);
        Vavg_3(lll) = deltaZ3(lll)./deltaT3(lll);

    end

    %display(['Simulation time is ' num2str(simulationTime/60/60) ' hours'])

end

if plotLegend==1 && plotBubbleArea==1
   
    figure(2)
    legend(cellLegend,'Location','Best')
    
    figure(3)
    legend(cellLegend,'Location','Best')
    
end

if plotLegend==1 && plotBubbleXcm==1
    
    figure(4)
    legend(cellLegend,'Location','Best')
    
    figure(5)
    legend(cellLegend,'Location','Best')
    
end

if plotLegend==1 && plotConeVel==1
    
    figure(6)
    legend(cellLegend,'Location','Best')
    
    figure(7)
    legend(cellLegend,'Location','Best')
    
    figure(40)
    legend(cellLegend,'Location','Best')
    
end

%plot final cone displacement
% figure
% loglog(param1,abs(finalPos),'o-')
% xyLabelTex(param1name,'\Delta z_c')
% %xlabel(param1name)
% %ylabel('T_{exit}')
% grid on
% title('Final cone position')

%plot coefficient of migration phase
figure
loglog(param1,plotMaxExcessArea,'ok-')
xyLabelTex(param1name,'||\Delta A||_\infty')
grid on

%plot coefficient of migration phase
figure
semilogx(param1,coeffPL,'ok-')
xyLabelTex(param1name,'\kappa')
grid on

%plot exponent of migration phase
figure
semilogx(param1,exponentPL,'ok-')
xyLabelTex(param1name,'\alpha')
grid on

%plot avergae velocty based on scaling
figure
vAvgScaling = coeffPL./(exponentPL+1).*deltaTmigration.^exponentPL;
loglog(param1,vAvgScaling,'ok-')
xyLabelTex(param1name,'\bar{V}_{scaling}')
grid on

%plot time to exit the cone
figure
loglog(param1,Texit,'o-')
xyLabelTex(param1name,'T_{\rm{exit}}')
hold on
loglog(param1,100*param1.^(-1/2),'-k')
legend('simulations','Ca^{-1/2}','Location','Best')
%xlabel(param1name)
%ylabel('T_{exit}')
grid on
title('Time to exit the cone')

figure
loglog(param1,1./Texit,'o-')
hold on
loglog(param1,0.05*param1.^0.66,'k-')
xyLabelTex(param1name,'f_{\rm{exit}}')
legend('simulations','Ca^{2/3}','Location','Best')
%xlabel(param1name)
%ylabel('T_{exit}')
grid on
title('Ejection frequency')

%plot final cone displacement
% figure
% plot(param1,TbubblePartiallyOut,'o-')
% hold on
% plot(param1,Texit-TbubblePartiallyOut,'o-')
% plot(param1,Texit,'o-')
% xyLabelTex(param1name,'\Delta t')
% %xlabel(param1name)
% %ylabel('T_{exit}')
% legend('\Delta t_1','\Delta t_2','\Delta t','Location','Best')
% grid on
% %title('Final cone position')
% 
% %plot ejection frquency
% figure
% plot(param1,1./TbubblePartiallyOut,'o-')
% hold on
% plot(param1,1./(Texit-TbubblePartiallyOut),'o-')
% plot(param1,1./(Texit),'o-')
% xyLabelTex(param1name,'f')
% %xlabel(param1name)
% %ylabel('T_{exit}')
% legend('f_1','f_2','f','Location','Best')
% grid on
% 
% %plot final cone displacement
% figure
% plot(param1,abs(ZCbubblePartiallyOut),'o-')
% hold on
% plot(param1,abs(finalPos-ZCbubblePartiallyOut),'o-')
% plot(param1,abs(finalPos),'o-')
% xyLabelTex(param1name,'|\Delta z_c|')
% %xlabel(param1name)
% %ylabel('T_{exit}')
% legend('\Delta z_{c}^1','\Delta z_{c}^2','\Delta z_{c}','Location','Best')
% grid on
% %title('Final cone position')
% 
% %plot average velocity
% figure
% plot(param1,abs(ZCbubblePartiallyOut./TbubblePartiallyOut),'o-')
% hold on
% plot(param1,abs(finalPos-ZCbubblePartiallyOut)./(Texit-TbubblePartiallyOut),'o-')
% plot(param1,abs(finalPos./Texit),'o-')
% xyLabelTex(param1name,'V_{\rm{avg}}')
% %xlabel(param1name)
% %ylabel('T_{exit}')
% legend('V_{\rm{avg}}^1','V_{\rm{avg}}^2','V_{\rm{avg}}','Location','Best')
% grid on
% title('average velocity')

%plot phases duration
figure
loglog(param1,deltaT1,'o-')
grid on
hold on
loglog(param1,deltaT2,'o-')
loglog(param1,deltaT3,'o-')
xyLabelTex(param1name,'\Delta T')
legend('\Delta T_{I}','\Delta T_{II}','\Delta T_{III}','Location','Best')
title('Phases duration')

%plot displacement in aech phase
figure
loglog(param1,deltaZ1,'o-')
grid on
hold on
loglog(param1,deltaZ2,'o-')
loglog(param1,deltaZ3,'o-')
xyLabelTex(param1name,'|\Delta Z|')
legend('|\Delta Z_I|','|\Delta Z_{II}|','|\Delta Z_{III}|','Location','Best')
title('Displacement in each phase')

%plot displacement in aech phase
figure
loglog(param1,Vavg_1,'o-')
grid on
hold on
loglog(param1,Vavg_2,'o-')
loglog(param1,Vavg_3,'o-')
xyLabelTex(param1name,'V^{avg}')
legend('V^{avg}_I','V^{avg}_{II}','V^{avg}_{III}','Location','Best')
title('Average vel in each phase')

%plot final cone displacement
figure
loglog(param1,abs(finalPos),'o-')
xyLabelTex(param1name,'z_c(\Delta T)')
hold on
loglog(param1,param1.^(1/4),'k-')
legend('simulations','Ca^{1/4}','Location','Best')
%xlabel(param1name)
%ylabel('T_{exit}')
grid on
title('Final cone position')

if numel(CaUP)>1 && numel(thetaUp)==1

    %plot avergae velocity
    figure
    loglog(param1,abs(vAverage),'o-')
    xyLabelTex(param1name,'V_{\rm{avg}}')
    hold on
    xxx = logspace(-3,-1,100);
    plot(xxx,0.01*xxx.^0.75,'k-')
    %xlabel(param1name)
    grid on
    %ylabel('V_{avg}')
    title('Average velocity')
    hold on
    loglog(param1,vAvgScaling,'or-')
    legend('From simulation','~Ca^{0.75}','k/(\alpha+1) (t_2-t_1)^\alpha','Interpreter','Latex')
    
    figure
    loglog(param1,abs(finalPos)./(param1'.*Texit),'o-')
    xyLabelTex(param1name,'\Delta z_c/(\rm{Ca} \Delta T)')
    %xlabel(param1name)
    %ylabel('T_{exit}')
    grid on
    title('Displacement efficiency')
    
    hold on
    xxx = logspace(-3,-1,100);
    plot(xxx,0.01*xxx.^-0.25,'k-')
    legend('From simulation','~Ca^{-0.25}','Location','Best')

    % fit power law
%     hold on
%     fitPowerLaw = fit(param1',abs(vAverage),'power1');
%     xxx = logspace(-5,-2,100);
%     yyy = fitPowerLaw.a*xxx.^fitPowerLaw.b;
%     loglog(xxx,yyy,'--k')
%     legend('From simulation',['Extrapolated ~' num2str(fitPowerLaw.b)],'Location','Best')

else
    
    %plot avergae velocity
    figure(20)
    if kkk==2
        hold on
    end
    plot(param1,abs(vAverage),'o-')
    xyLabelTex(param1name,'V_{\rm{avg}}')
    %xlabel(param1name)
    grid on
    %ylabel('V_{avg}')
    title('Average velocity')
    %legendParam2()
    
end

if plotDimRadius==1

    figure
    loglog(param1/param1(1)*Lscale,abs(vAverage)*Vscale,'o-')
    xyLabelTex('R','V_{\rm{avg}} [\rm{m}/\rm{s}]')
    %xlabel(param1name)
    grid on
    %ylabel('V_{avg}')
    title('Like increasing R')

end

if numel(param2)>1

    [maxAvgVel(kkk),ind] = max(abs(vAverage));
    thetaMaxAvgVel(kkk) = thetaUp(ind);
    
    Vavg1 = abs(ZCbubblePartiallyOut)./(TbubblePartiallyOut);
    Vavg2 = abs(finalPos-ZCbubblePartiallyOut)./(Texit-TbubblePartiallyOut);
    
    [maxAvgVel1(kkk),ind1] = max(abs(Vavg1));
    thetaMaxAvgVel1(kkk) = thetaUp(ind1);
    
    [maxAvgVel2(kkk),ind2] = max(abs(Vavg2));
    thetaMaxAvgVel2(kkk) = thetaUp(ind2);

end

end

if numel(param2)>1

figure
loglog(param2,maxAvgVel,'o-k')
xyLabelTex(param2name,' \rm{max} ( V_{\rm{avg}}) ')
grid on
title('Maximum average velocity')

figure
plot(param2,thetaMaxAvgVel,'o-k')
xyLabelTex(param2name,' \theta^{\rm{opt}} ')
grid on
title('Optimal theta')

figure
plot(param2,thetaMaxAvgVel1,'o-k')
xyLabelTex(param2name,' \theta^{\rm{opt}}_1 ')
grid on
title('Optimal theta 1')

figure
plot(param2,thetaMaxAvgVel2,'o-k')
xyLabelTex(param2name,' \theta^{\rm{opt}}_2 ')
grid on
title('Optimal theta 2')

end













