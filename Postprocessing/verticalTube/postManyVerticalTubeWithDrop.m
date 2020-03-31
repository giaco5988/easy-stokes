%post processing of a conical motor with a defromable bubble

clear variables
close all

%parameters
res = 0;
%BoUP = [2 1.8 1.6 1.4 1.2 1 0.95 0.9];
%BoUP = [2 1.5 1.3 1 0.8 0.5 0.1 0.05 0.04 0.03 0.02 0.01];
BoUP = 0.2:0.2:2;
%BoUP = [0.2:0.2:1 1.6 1.8 2];
%BoUP = 1.2;
lambdaUP = 0;
dtUp = 0.02;
%TendUPpost = [5e2 1e3 1e3 1e3 1.1e3 2e3 2e3 2e3];
TendUPpost = 2001;
%TendUPpost = [2000 2000 2000 2000 2000 500 500 500];
alphaUPpost = 1.4;
%alphaUPpost = [0.6 0.8:0.05:1.4];
ODE = 0;
%Lup = 20:20:100;
Lup = 20;
%nElemUP = [20 40 60 80 100]+5*Lup+10;
nElemUP = 300;
%nElemUP = round(alphaUPpost*50)+10*Lup+20;
%x0Upload = -Lup/2+5;
x0Upload = 0;
repUPpost = 5;

param1 = BoUP;  param1name = 'Bo';
%param1 = Lup;  param1name = '\xi';
%param1 = x0Upload;  param1name = 'x_0';
%param1 = nElemUP;  param1name = 'N';
%param1 = alphaUPpost;  param1name = '\alpha';

color = get(gca,'ColorOrder');

%options
plotLegend = 1;
plotRes = 0;
plotDropVel = 1;    notLastState = [1 1 1 1 1 10 10 20 50 20];   %notLastState = [30 ones(1,100)];
plotFilm = 1;
plotBubbleArea = 0;
plotBubbleVolume = 0;
plotBubbleXcm = 0;
plotLastCurvature = 1;
plotVelField = 0;   Tplot = 24;    substitutePoint = 1;    coeffSub = 0.5;
plotVelFieldSlice = 0;  xSlice = 2; nSlice = 100;   Ytop = 1.4;
cutX = -3;   shift = 0;  cutY = 1;   nx = 10;    ny = 10;
plotVelDecay = 0;   nDecay = 5;
compareWithLudo = 1;    pathLudo = '/Users/Giacomo/Documents/research_notes/verticalPipe/compareDataMeAndLudo/';

manyLastFilm = zeros(numel(param1),1);
manyLastDropVel = zeros(numel(param1),1);
curvTipRight = zeros(numel(param1),1);
curvTipLeft = zeros(numel(param1),1);
cellLegend = cell(numel(param1),1);
for kkk = 1:numel(param1)

BoUPhere = BoUP(kkk);
LupHere = Lup; nElemUPhere = nElemUP; x0UploadHere = x0Upload;
dtUpHere = dtUp;
TendUPpostHere = TendUPpost;
alphaUPhere = alphaUPpost;

if plotLegend==1
     cellLegend{kkk} = [param1name '=' num2str(param1(kkk))];
end
    
%filename
filename = ['verticalTube_ODE=' num2str(ODE) '_rep=' num2str(repUPpost) '_x0=' num2str(x0UploadHere) '_alpha=' num2str(alphaUPhere) '_Bond=' num2str(BoUPhere) '_lambda=' num2str(lambdaUP) '_L=' num2str(LupHere) '_nStart=' num2str(nElemUPhere) '_Tend=' num2str(TendUPpostHere) '_dt=' num2str(dtUpHere) '.mat'];

%source
if res==1
   source = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==0
   source = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==2
   source = '~/Documents/MATLAB/test/verticalPipeDrop/resTest/';
end

BEM = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/verticalTube';
here = pwd;

saveDest = '/Users/Giacomo/Documents/research_notes/APS_2017/movies/frames/';
namePlot = 'plotConeVelocityDeformable';

%load file
load([source filename])
PARAM.kernelFreeSpace = 1;  PARAM.posWall = [];
PARAM.ellipseShape = [0 0 0 1];
PARAM.D = [0 0 0 0];

%number
if T(end)==0
    loopUp = find(T==0,2);
    loopUp = loopUp(2)-1;
else
    loopUp = numel(T);
end
%initialize
manyRes = zeros(numel(T),1);
hFilmPost = zeros(numel(T),1);
minDL = zeros(numel(T),1);
dropVel = zeros(numel(T),1);
frontVel = zeros(numel(T),1);
Area = zeros(numel(T),1);
Volume = zeros(numel(T),1);
xcm = zeros(numel(T),1);
maxBubble = zeros(numel(T),1);
tParametricPost = {linspace(0,1,10) linspace(0,1,100) linspace(0,1,10) nan};
for i = 1:loopUp
    
   display([num2str(i) ' of ' num2str(loopUp)])
    
   %current shape and motor position
   Yhere = Y{i};
   if ODE==0
    Vhere = V{i};
   end
   xBubble = Yhere(1:2:end-1);
   yBubble = Yhere(2:2:end);
   
   %current location   TEMPORARY!!!
   PARAM_now = PARAM;
   
   %compute curent shape
   [x,y] = buildGeometryPanelsParametric(tParametricPost,PARAM_now);
   x{4} = xBubble';
   y{4} = yBubble';
   
   %residuals
   if plotRes==1 || plotDropVel==1
    if ODE~=0
        cd(BEM)
        Vhere = computeVelocityRising(T(i),Yhere,tParametricBase,PARAM_now);
        cd(here)
    end
    [nxBubble,nyBubble] = computeNormalVector(xBubble',yBubble',PARAM.orderVariableStokes(4),PARAM.orderGeometryStokes(4),PARAM.SPlinesType(4));
    dropVel(i) = DropVelocityAxis(xBubble',yBubble',Vhere(1:2:end-1).*nxBubble'+Vhere(2:2:end).*nyBubble');
    frontVel(i) = Vhere(1);
    res = NormalVelocityDropFrame(Yhere(1:2:end-1)',Yhere(2:2:end)',Vhere(1:2:end-1),Vhere(2:2:end));
    manyRes(i) = norm(res,Inf);
   end
   
   %bubble surface area
   Area(i) = surf_gauss_vect(xBubble',yBubble');
   
   %bubble volume
   Volume(i) = axis_int_gauss_vect(xBubble',yBubble');
   
   %bubble center of mass
   xcm(i) = center_mass(xBubble',yBubble');
   maxBubble(i) = max(xBubble);
   %xcm(i) = centerOfMassBlockAxis(x,y,2,PARAM);
   
   %compute distance wall-drop
   %hFilmPost(i) = min(distWallDrop(x{4},y{4},x{2},y{2}));
   hFilmPost(i) = min(1-y{4});
   dx = diff(xBubble);
   dy = diff(yBubble);
   dl = sqrt(dx.^2+dy.^2);
   minDL(i) = min(dl);
    
end

manyLastFilm(kkk) = hFilmPost(i-notLastState(kkk));
Tlast(kkk) = T(i-notLastState(kkk));

if plotDropVel==1
    
    %kkk
    manyLastDropVel(kkk) = dropVel(i-notLastState(kkk));
    
end

%derivative
D1 = finiteDifference1D(i,[2 0],1);

%area for a sphere
Rarea = nthroot(3/4/pi*Volume,3);
AreaSphere = 4*pi*Rarea;

if plotLastCurvature==1
    
   %curvature
   PARAM_now.n(4) = numel(x{4})-1;
   %[nx,ny] = computeNormalVector(x{4},y{4},PARAM_now.orderVariableStokes(4),PARAM_now.orderGeometryStokes(4),PARAM_now.SPlinesType(4));
   [k1,k2] = computeCurvatureSplines(x{4},y{4},PARAM_now.orderVariableStokes(4));
   k = k1+k2;
    
   figure(2)
   dx = diff(x{4}); dy = diff(y{4}); dl = sqrt(dx.^2+dy.^2);
   %lArc = [0 cumsum(dl)];
   %plot(lArc,k,'-x')
   if kkk>1
        hold on
    end
   plot(x{4},k,'-')
   xlabel('z')
   ylabel('K')
   title('Curvature')
   grid on
   curvTipRight(kkk) = k(1);
   curvTipLeft(kkk) = k(end);
    
end

%plot bubble volume and its variation
if plotBubbleVolume==1
    
    figure(3)
    %subplot(2,1,1)
    %hold on
    Vstart = 4/3*pi*alphaUPhere^3;
    err = abs(Volume(1:i)-Vstart)/Vstart;
    semilogy(T(1:i),err)
    grid on
    xlabel('t')
    ylabel('errV')
    title('Error on Volume')
    drawnow
      
end

%plot bubble surface area
if plotBubbleArea==1
    
    figure
    %subplot(2,1,1)
    plot(T(1:i),Area(1:i)-AreaSphere(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('A')
    title('Excess surface area')
    drawnow
    
    dAdt = (D1*Area(1:i))./(D1*T(1:i));
    dAdtSphere = (D1*AreaSphere(1:i))./(D1*T(1:i));
    
    %subplot(2,1,2)
    figure
    plot(T(1:i),dAdt(1:i)-dAdtSphere(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('dA/dt')
    title('Surface area variation')
    drawnow
      
end

%plot bubble center of mass
if plotBubbleXcm==1
    
    figure
    %subplot(2,1,1)
    plot(T(1:i),xcm(1:i),'-k')
    grid on
    xlabel('t')
    ylabel('X_{cm}')
    title('Bubble position')
    drawnow
      
end

%plot film thicjness
if plotFilm==1
    
    figure(20)
    ind = find(xcm==0,1,'first');
    if isempty(ind)==1
        plot(xcm(end)+Lup/2,hFilmPost(end),'.k','MarkerSize',30)
    else
        plot(xcm(ind)+Lup/2,hFilmPost(ind),'.k','MarkerSize',30)
    end
    hold on
    
    figure(7)
    if kkk>1
        hold on
    end
    %semilogy(T(1:i),hFilmPost(1:i))
    plot(T(1:i),hFilmPost(1:i))
    grid on
    xlabel('t')
    ylabel('h')
    title('film thickness')
    hold on
    plot(Tlast(kkk),manyLastFilm(kkk),'.k','MarkerSize',30)
    drawnow
    
    figure(10)
    if kkk>1
        hold on
    end
    %semilogy(T(1:i),hFilmPost(1:i))
    plot(xcm(1:i)+LupHere/2,hFilmPost(1:i))
    grid on
    xlabel('x_{cm}')
    ylabel('h')
    title('film thickness')
    drawnow
    
    %flim thickness variation
    varFilm(kkk) = (hFilmPost(i)-hFilmPost(round(i/4*3)))/(T(i)-T(round(i/4*3)));
       
end

%plot bubble velocity
if plotRes==1
    
    figure(8)
    if kkk>1
        hold on
    end
    %hold on
    semilogy(T(1:i),manyRes(1:i))
    %plot(T(1:i),manyRes(1:i))
    grid on
    xlabel('t')
    ylabel('R')
    title('Residuals')
    drawnow
       
end

if plotDropVel==1
    
    figure(9)
    if kkk>1
        hold on
    end
    %semilogy(T(1:i),dropVel(1:i))
    plot(T(1:i),dropVel(1:i))
    grid on
    xyLabelTex('t','u_d')
    title('Droplet velocity')
    hold on
    plot(Tlast(kkk),manyLastDropVel(kkk),'.k','MarkerSize',30)
    drawnow
    
    figure(11)
    if kkk>1
        hold on
    end
    %semilogy(T(1:i),dropVel(1:i))
    plot(maxBubble(1:i),frontVel(1:i))
    grid on
    xyLabelTex('z_{front}','u_{front}')
    title('Front velocity')
    drawnow
       
end

disp(['Simulation time is ' num2str(simulationTime/60/60) ' hours'])

end

%plot legend
if plotLegend==1

    if plotDropVel==1
        figure(9)
        legend(cellLegend,'Location','Best')
    end
    
    if plotRes==1
        figure(8)
        legend(cellLegend,'Location','Best')
    end
    
    if plotFilm==1
        figure(7)
        legend(cellLegend,'Location','Best')
    end
    
    if plotLastCurvature==1
        figure(2)
        legend(cellLegend,'Location','Best')
    end

end

%plot film thicness
figure
%loglog(param1,manyLastFilm,'x','MarkerSize',20)
semilogy(param1,manyLastFilm,'x','MarkerSize',20)
%xlabel(param1name)
grid on
%ylabel('h')
xyLabelTex(param1name,'h')
title('Film thickness')

if compareWithLudo==1
    
    load([pathLudo 'Ca_and_h_flat_h_rear_and_r_cap_vs_Bo.mat'])
    hold on
    plot(Bo,h_minimal_rear./r_cap,'ok')
    legend('BEM (Giac)','Comsol (Ludo)','Location','Best')
    
end

%plot last drop vel
figure
%plot(param1,manyLastDropVel,'x','MarkerSize',20)
semilogy(param1,manyLastDropVel,'x','MarkerSize',20)
%loglog(param1,manyLastDropVel,'x','MarkerSize',20)
%xlabel(param1name)
grid on
%ylabel('u_d')
xyLabelTex(param1name,'u_d')
title('Droplet velocity')

if compareWithLudo==1
    
    load([pathLudo 'Ca_and_h_flat_h_rear_and_r_cap_vs_Bo.mat'])
    hold on
    plot(Bo,Ca,'ok')
    legend('BEM (Giac)','Comsol (Ludo)','Location','Best')
    
end

if numel(BoUP)>1
    %plot last drop vel
    figure
    loglog(param1-0.82,manyLastDropVel,'x','MarkerSize',20)
    xlabel(param1name)
    grid on
    ylabel('u_d')
    xyLabelTex('Bo-Bo_c','u_d')
    title('Droplet velocity')
    
    if compareWithLudo==1
    
    load([pathLudo 'Ca_and_h_flat_h_rear_and_r_cap_vs_Bo.mat'])
    hold on
    plot(Bo-0.82,Ca,'ok')
    legend('BEM (Giac)','Comsol (Ludo)','Location','Best')
    
    end
    
end

%plot curv tip
figure
plot(param1,curvTipRight,'-')
hold on
plot(param1,curvTipLeft,'-')
xlabel(param1name)
grid on
ylabel('k')
legend('k right','k left','Location','Best')
title('Curvature of the caps')

figure
plot(param1,varFilm,'.-','MarkerSize',30)
xlabel(param1name)
grid on
ylabel('dh/dt')
%legend('k right','k left','Location','Best')
%title('Film variation in time')














