%post processing of test time stepping spectral

close all
clear variables

%% Add libraries paths
REPOSITORY_NAME = '~/Documents/MATLAB/';    % path to the repository
addpath('../utils_one_drop_spectral')
add_paths_tutorials_drop_spectral(REPOSITORY_NAME);

%% Parameters of simulation to upload
Ca = 0.1;       % capillary number
visc = 5;       % viscosity ratio
n = 25;         % number of modes
Tend = 500;     % time of simulation
ODE = 2;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1;         % 1 is extensional flow, 2 is rising droplet
volCorr = 0;    % if volume correction has been used
maxDT = 1e-2;   % max time step
legendre = 1;   % if 1, use Legendre, if 0, Chebyshev
D = 0.0;        % deformation parameter of the initial shape

%% Plotting options
plotShape = 1;  threeD = 1;     % if 1, plot droplet shape, and optionally use 2d rendering
plotRes = 1;                    % if 1, plot residuals
Tbreak = 2*Tend;                % time to stop the post-processing

%% Upload data
dir = '../tutorial_results/';
name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(D) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(Tend) '.mat'];
load([dir name])

% %initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%% Varibale initialization
ite = round(numel(T));
ForceFree = zeros(ite,1);
normElon = zeros(ite,1);
res = zeros(ite,1);
V = zeros(ite,1);
Area = zeros(ite,1);
Evisc = zeros(ite,1);
Dellipse = zeros(ite,1);
xcm = zeros(ite,1);
ExcessArea = zeros(ite,1);
k1_east = zeros(ite,1);
k1_north = zeros(ite,1);
k2_east = zeros(ite,1);
k2_north = zeros(ite,1);
f1 = zeros(ite,1);
f2 = zeros(ite,1);
f3 = zeros(ite,1);
V0 = 4/3*pi;

%loops
plotCount = 1;
shift = 0;  islast = 0;
for i = 1:ite
    
   k = i;
    
   display([num2str(i) ' of ' num2str(ite)])
    
   %get currents modes
   xyMode = Y(k,:)';
   xMode = Y(k,1:2:end-1)';
   yMode = Y(k,2:2:end)';
   
   %interpolation
    if PARAM.legendre==2 || PARAM.legendre==1

          xInterp = @(t) interpLegendreZeroOne(t,xMode);
          yInterp = @(t) interpLegendreZeroOne(t,yMode);

    elseif PARAM.legendre==0

          xInterp = chebcoeffs2chebvals(xMode);
          yInterp = chebcoeffs2chebvals(yMode);
          
          xInterp = chebfun(xInterp,[0 1]);
          yInterp = chebfun(yInterp,[0 1]);
          
    end
   
   [x,y] = fromModesToGrid(xMode,yMode,PARAM);
   
   %compute elongation norm
   RforNorm = sqrt(x.^2+y.^2);
   xcm(i) = CenterMassCurvAxisSpectral(x,y,PARAM);
   if Ca>=0
        normElon(i) = abs(RforNorm(1)-xcm(i));
   else
       [~,ind] = max(abs(y));
        normElon(i) = y(ind);
   end
   
   %compute volume abd surface area
   V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);
   Area(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM);
   R0 = nthroot(V(i)/4/pi*3,3);
   ExcessArea(i) = Area(i) - 4*pi*R0^2;
   
   %defromation parameter
   L = xInterp(0)-xInterp(1);
   B = 2*yInterp(0.5);
   Dellipse(i) = (L-B)/(L+B);
   
   %compute curvature
   xp = PARAM.D1*x;    xpp = PARAM.D2*x;
   yp = PARAM.D1*y;    ypp = PARAM.D2*y;
   h = (xp.^2+yp.^2).^(0.5);
   nx = yp./h;
   ny = -xp./h;
   K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
   if PARAM.legendre==1
        K2 = ny./y;
   elseif PARAM.legendre==0||PARAM.legendre==2
        K2 = ny./y;
        K2([1 end]) = K1([1 end]);
   else
       error('not implemented')
   end
   
   if plotShape==1
       
      xMid = (max(x)+min(x))/2;
      
      figure(1)
      if threeD==1
        
          solidOfRevolution(y',x',[0 2*pi],[0 0.8 1],0.8)
          axis off
          axis equal
          axis([-5 5 -1 1 -1 1])
        
      else
          if BC==2
                plot(y,-x,'k-')
                hold on
                grid on
                plot(-y,-x,'k')
                title('Rising Droplet')
                axis([-2 2 -2 2])
          elseif BC==1 || BC==3
                plot(x,y,'k-')
                hold on
                grid on
                plot(x,-y,'k')
                axis([-4 4 -2 2])
                title('Extensional')
          end
          xlabel('r')
          ylabel('x')
          axis equal
          %axis([-2 2 -2-xMid 2-xMid])
          hold off
          %title('Rising Droplet')
      end
      
      if zoom==1
          
         axis([-0.1 0.1 -1.2 -1.0]) 
          
      end
      
      drawnow
      
   end
   
   if plotRes==1
           
           if PARAM.legendre==1||PARAM.legendre==2
               [~,res(i)] = tutorial_dropLegendreCurvilinearModes(0,xyMode,PARAM);
           elseif PARAM.legendre==0
               [~,res(i)] = tutorial_dropExtensChebfunCurvilinearModes(0,xyMode,PARAM);
           end
       
   end
   
   %stop when there are no more data
   if i~=ite
       if T(k+1)>Tbreak||(T(k+1)==0&&k>1)
           ForceFree = ForceFree(1:i);
           V = V(1:i);
           t = T(1:i);
           res = res(1:i);
           normElon = normElon(1:i);
           disp('Break')
           break;
       end
   else
       t = T;
   end
    
end

figure
plot(t,normElon,'k-')
grid on
title('Elongation norm')
xlabel('t')
ylabel('||r(0)-1||_{\infty}')

figure
loglog(abs(xMode(2:2:end)))
hold on
loglog(abs(yMode(1:2:end)))
grid on
hold off
xlabel('n')
ylabel('c')
legend('C_x','C_y')
title('Modes coefficients')

errV = (V-V0)/V0;
figure
plot(t,errV)
xlabel('t')
ylabel('err_V')
grid on
title('Error On Volume')

figure
plot(t(1:i),ExcessArea(1:i))
xlabel('t')
ylabel('\Delta A')
title('Excess area')
grid on

figure
semilogy(t,res)
hold on
xyLabelTex('t','|\bf u \cdot n|')
title('residuals')
grid on

%compute curvature
xp = PARAM.D1*x;    xpp = PARAM.D2*x;
yp = PARAM.D1*y;    ypp = PARAM.D2*y;
h = (xp.^2+yp.^2).^(0.5);
nx = yp./h;
ny = -xp./h;
K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
if PARAM.legendre==1
    K2 = ny./y;
elseif PARAM.legendre==0||PARAM.legendre==2
    K2 = ny./y;
    K2([1 end]) = K1([1 end]);
end

%subplot(3,2,6)
figure
plot(x,K1)
hold on
plot(x,K2)
plot(x,K1+K2,'k')
xlabel('x')
ylabel('k')
legend('k_1','k_2','k','Location','Best')
title('curvature')
grid on

%% Print simulation time
disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











