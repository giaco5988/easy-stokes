%post processing of test time stepping spectral

close all
clear variables

%% ADD LIBRARIES
REPOSITORY_NAME = '~/Documents/MATLAB/';    % path to the repository
addpath('../utils_one_drop_spectral')
add_paths_tutorials_drop_spectral(REPOSITORY_NAME);

%% PATH TO RESULTS
dir = '../tutorial_results/';

%% UPLOAD PARAMETERS
Ca = 0.1;               % Capillary number
visc = 5;               % viscosity ratio
n = 50;                 % space discretization
ODE = 2;                % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1;                 % 1 is extensional flow, 2 is rising droplet
volCorr = 0;            % volume correction flag
maxDT = 0.01;           % maximum time step if adaptive, otherwise simply time step
legendre = 1;           % if 1 use Legendre, if 0 use Chebyshev
edgeLoop = 20;          % number of edge tracking iteration
edgeDelta = 1e-2;       % bisect trajectory when they are this far

%% PLOTTING OPTION
colIN = '--'; colOUT = '-';     % line for stable/unstable trajectories
plotShape = 1;   mesh = 0;      % if 1, plot shape
plotNorm = 1;                   % if 1, plot droplet elongation
plotRes = 1;                    % if 1, plot residuals
plotVol = 0;                    % if 1, plot volume error
plotPhaseSpacePolar = 0;    symmetric = 0;  modes = 7;  % if 1, plot trajecotries in phase space
convergeShape  = 0.02;          % set criterion to determine convergence

%% UPLOAD DATA
name = ['edgeTrackingDrop_edgeLoop=' num2str(edgeLoop) '_deltaEdge=' num2str(edgeDelta) '_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_maxDT=' num2str(maxDT) '_VolCorr=' num2str(volCorr) '.mat'];
load([dir name])

%% INITIALIZE CHECBFUN
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%get colot rgb
Tbreak = Inf;
step = 1;
color = get(gca,'ColorOrder');
close gcf
resEdge = 1;
for l = 1:edgeLoop
    
    T = Tedge{l};
    Y = Yedge{l};
    display([num2str(l) ' of ' num2str(numel(Tedge))])
    
    if isempty(T)
        warning('No data')
        break;
    end

    %number of iteration
    ite = round(numel(T)/step);

    %initialize
    ForceFree = zeros(ite,1);
    normElon = zeros(ite,1);
    normL2 = zeros(ite,1);
    res = zeros(ite,1);
    V = zeros(ite,1);
    DDD = zeros(ite,1);
    A = zeros(ite,1);
    V0 = 4/3*pi;
    f1 = zeros(ite,1);
    f2 = zeros(ite,1);
    f3 = zeros(ite,1);

    %loops
    countShift = 1;
    countShape = 1;
    for i = 1:ite

       k = i*step;

       %get currents modes
       xyMode = Y(k,:)';
       xMode = Y(k,1:2:end-1)';
       yMode = Y(k,2:2:end)';

       if PARAM.legendre==1||PARAM.legendre==2

           x = LegendreBuildXY(xMode,PARAM.PPP);
           y = LegendreBuildXY(yMode,PARAM.PPP);

       elseif PARAM.legendre==0

           x = chebcoeffs2chebvals(xMode);
           y = chebcoeffs2chebvals(yMode);

       end
       
       if plotPhaseSpacePolar==1
            
            %radius modes
            f = LegendreSerie(x,y,modes,symmetric,PARAM);
            f1(i) = f(3);
            f2(i) = f(5);
            f3(i) = f(7);
        
       end

       %compute volume
       V(i) = VolumeCurvilinearAxisSpectral(x,y,PARAM);
       xcm = CenterMassCurvAxisSpectral(x,y,PARAM);
       R0 = nthroot(3/4*V(i)/pi,3);
       
       %deformation parameter
       L = max(x)-min(x);
       B = 2*y(round(numel(y)/2));
       DDD(k) = (L-B)/(L+B);
       
       %compute elongation norm
       RforNorm = sqrt(x.^2+y.^2);
       normElon(i) = x(1)-xcm;
       %normElon(i) = DDD(k);
       normL2(i) = sqrt(PARAM.WG'*(RforNorm-1).^2);
       
       %compute area
       A(i) = surfaceCurvilinearAxisSpectral(x,y,PARAM);

       %compute force acting on droplet
       ForceFree(i) = ForceOnDropSpectralCurvilinear(x,y,PARAM.Ca,PARAM);

       t = T;

       if plotShape==1

          xMid = (max(x)+min(x))/2;

          figure(1) 
          if BC==2
                plot(y,-x,'k-')
                hold on
                grid on
                plot(-y,-x,'k')
                title('Rising Droplet')
                axis([-2 2 -2-xcm 2-xcm])
                if mesh==1
                  plot(y,-x,'xr')
                end
                xlabel('r')
          ylabel('x')
          elseif BC==1
                plot(x,y,'k-')
                hold on
                grid on
                plot(x,-y,'k')
                axis([-3 3 -3 3])
                title(['Edge shape Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
                if mesh==1
                  plot(x,-y,'xr')
                end
                xlabel('z')
                ylabel('r')
          end
          axis equal
          %title('Rising Droplet')
          hold off
          drawnow
          
       end
       
       if plotRes==1

               if PARAM.legendre==1
                   [~,res(i)] = tutorial_dropLegendreCurvilinearModes(t,xyMode,PARAM);
               elseif PARAM.legendre==0
                   [~,res(i)] = tutorial_dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
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
               normL2 = normL2(1:i);
               disp('Break')
               break;
           end
       end

    end
    
    finalInd = numel(t);
    
    %figure out if it is in or out the bassin of attraction
    L = max(x)-min(x);
    B = 2*y(round(numel(y)/2));
    D = (L-B)/(L+B);
    
    if PARAM.BC==1
    
        if PARAM.visc==1
                    load('./steadyState/CaExt')
                    load('./steadyState/DExt')
        elseif PARAM.visc==0
                    load('./steadyState/CaExt0')
                    load('./steadyState/DExt0')
        elseif PARAM.visc==0.02
                    load('./steadyState/CaExt002')
                    load('./steadyState/DExt002')
        elseif PARAM.visc==0.05
                    load('./steadyState/CaExt005')
                    load('./steadyState/DExt005')
        elseif PARAM.visc==0.1
                    load('./steadyState/CaExt01')
                    load('./steadyState/DExt01')
        elseif PARAM.visc==0.5
                    load('./steadyState/CaExt05')
                    load('./steadyState/DExt05')
        elseif PARAM.visc==5
                    load('./steadyState/CaExt5')
                    load('./steadyState/DExt5')
        elseif PARAM.visc==10
                    load('./steadyState/CaExt10')
                    load('./steadyState/DExt10')
        elseif PARAM.visc==0.01
                    load('./steadyState/CaExt001')
                    load('./steadyState/DExt001')
        end
        [~,ind] = min(abs(PARAM.Ca-manyCa));
        Dca = manyD(ind);
        OutIn  = abs(D-Dca)/Dca > convergeShape;
    
    elseif PARAM.BC==2

        OutIn  = abs(D) > convergeShape;

    end
    
    if OutIn==1
        col = colOUT;
    elseif OutIn==0
        col = colIN;
    end
    
    if plotVol==1
        errV = (V-V0)/V0;
        figure(5)
        semilogy(t,errV)
        hold on
        xlabel('t')
        ylabel('err_V')
        grid on
        title('Error On Volume')
        drawnow
    end
    
    if plotRes==1
        
        figure(20)
        semilogy(t,res,col)
        hold on
        xlabel('t')
        ylabel('res')
        title(['residuals Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        grid on
        drawnow
        
    end

    if plotNorm==1
        figure(3)
            
        if BC==2
            semilogy(t(1:finalInd),normElon(1:finalInd),col)
            hold on
            grid on
            title('Rising Droplet')
            xlabel('t')
            ylabel('L/a')
        elseif BC==1
            semilogy(t(1:finalInd),normElon(1:finalInd),col)
            %semilogy(t,normElon,col)
            grid on
            hold on
            title(['Extensional flow Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
            xlabel('t')
            ylabel('L/2')
        end
        drawnow
        
    end
    
    %plot phase space
    if plotPhaseSpacePolar==1
        
        figure(11)
        
        if l>rangeLoop(1)
            hold on
        end
        plot(f1(1:finalInd),f2(1:finalInd),col)
        grid on
        
        xlabel('f_2')
        ylabel('f_4')
        title(['Phase space Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        
    end

end

%% print simulation time
disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])











