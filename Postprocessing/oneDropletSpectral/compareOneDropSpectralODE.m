%compare results of test time stepping spectral

close all
clear variables

%data
manyCa = 0.05*ones(4,1);
visc = [0 0.1 1 10];
n = [25 25 25 25 50 100 50 100 50 100 50 100 50 100 50 100];
manyTend = 100*ones(16,1);
manyODE = [2 2 2 2 1 1 1 1 6 6 6 6 1 1 1 1];                        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1*ones(16,1);                         % 1 is extensional flow, 2 is rising droplet
volCorr = 0*ones(1,4);
MANYmaxDT = 1*[2e-2 2e-2 2e-2 2e-2 1e-1 1e-1 1e-2 1e-2 1e-1 1e-1 1e-2 1e-2 1e-1 1e-1 1e-2 1e-2];
legendre = 1*ones(16,1);

%number of comparison
compare = numel(volCorr);

%option
plotShape = 0;
plotModes = 0;
plotRes = 1;
Tbreak = 200;
step = 1;

%directory
%dir = '~/Documents/MATLAB/droplet_simulations/server/';
dir = '~/Documents/MATLAB/droplet_simulations/results/';

%filename
name = ['DropSpectral_ODE=' num2str(manyODE(1)) '_Legendre=' num2str(legendre(1)) '_BC=' num2str(BC(1)) '_Ca=' num2str(manyCa(1)) '_visc=' num2str(visc(1)) '_n=' num2str(n(1)) '_maxDT=' num2str(MANYmaxDT(1)) '_volCorr=' num2str(volCorr(1)) '_Tend=' num2str(manyTend(1)) '.mat'];

%upload data
upload = [dir name];
load(upload)

%initialize chebfun
if PARAM.legendre==0
    x = chebfun('x',[0 1]);     y = chebfun('y',[0 1]);
end

%initialize
ForceFree = zeros(SaveHowMany,compare);
res = zeros(SaveHowMany,compare);
V = zeros(SaveHowMany,compare);
cellLegend = {zeros(1,compare)};
V0 = 4/3*pi;

%loop for comparison
for l = 1:compare
    
ODE = manyODE(l);   
 
if ODE==1
    ODEwrite = 'ODE45';
elseif ODE==6
    ODEwrite = 'ODE23t';
end
    
%print to screen
display([num2str(l) ' of ' num2str(compare)])
    
%filename
name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre(l)) '_BC=' num2str(BC(l)) '_Ca=' num2str(manyCa(l)) '_visc=' num2str(visc(l)) '_n=' num2str(n(l)) '_maxDT=' num2str(MANYmaxDT(l)) '_volCorr=' num2str(volCorr(l)) '_Tend=' num2str(manyTend(l)) '.mat'];

%upload data
upload = [dir name];
load(upload)

%number of iteration
ite = round(numel(T)/step);

%store element in cell for legend
cellLegend(l) = {[ODEwrite ' n=' num2str(n(l)) ' maxDT=' num2str(MANYmaxDT(l)) ' V=' num2str(volCorr(l)) ' t=' num2str(simulationTime)]};

%loops
for i = 1:ite
    
   k = i*step;
    
   %get currents modes
   if PARAM.volume==1
       if ODE==2
            xyMode = Y(:,k);
            xMode = xyMode([1 2:2:end-1]);
            yMode = xyMode(3:2:end);
       else
            xyMode = Y(k,:)';
            xMode = xyMode([1 2:2:end-1]);
            yMode = xyMode(3:2:end);
       end
       
       %compute modes first y mode from volume (x is set to zero)
       fVolume = @(rho) ModifyVolumeSpectralXYmodes(xMode,yMode,rho,V0,PARAM);
       options = optimoptions('fsolve','TolFun',1e-15,'TolX',1e-15,'Display','none');
       yFirst = fsolve(fVolume,1,options);
       
       yMode = [yFirst; yMode];
       
   elseif PARAM.volume==0
       
       if ODE==2
           xyMode = Y(:,k);
            xMode = Y(1:2:end-1,k);
            yMode = Y(2:2:end,k);
       else
           xyMode = Y(k,:)';
            xMode = Y(k,1:2:end-1)';
            yMode = Y(k,2:2:end)';
       end
   end
   
   if PARAM.legendre==1
       
       x = LegendreBuildXY(xMode,PARAM.PPP);
       y = LegendreBuildXY(yMode,PARAM.PPP);
       
   elseif PARAM.legendre==0
   
       x = chebcoeffs2chebvals(xMode);
       y = chebcoeffs2chebvals(yMode);
   
   end
   
   %compute volume
   V(i,l) = VolumeCurvilinearAxisSpectral(x,y,PARAM);
   
   %compute force acting on droplet
   ForceFree(i,l) = ForceOnDropSpectralCurvilinear(x,y,PARAM.Ca,PARAM);
   
   t = T;
   
   if plotShape==1
       
      xMid = (max(x)+min(x))/2;
       
      figure(1) 
      if BC==2
            plot(y,-x,'kx-')
            hold on
            grid on
            plot(-y,-x,'k')
            title('Rising Droplet')
      elseif BC==1
            plot(x,y,'kx-')
            hold on
            grid on
            plot(x,-y,'k')
            %axis([-2 2 -2 2])
            title('Extensional')
      end
      xlabel('r')
      ylabel('x')
      axis equal
      %axis([-2 2 -2-xMid 2-xMid])
      hold off
      %title('Rising Droplet')
      drawnow
      
   end
   
   if plotModes==1
        %plot modes
        figure(2)
        loglog(abs(xMode(2:2:end)))
        hold on
        loglog(abs(yMode(1:2:end)))
        grid on
        hold off
        xlabel('n')
        ylabel('c')
        legend('C_x','C_y')
        drawnow
   end
   
   if plotRes==1
       
       if PARAM.volume==1
           
           if PARAM.legendre==1
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [~,res(i,l)] = dropLegendreCurvilinearModesVolume(t,xyMode,V0,PARAM);
               cd(here)
           elseif PARAM.legendre==0
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [~,res(i,l)] = dropExtensChebfunCurvilinearModesVolume(t,xyMode,V0,PARAM);
               cd(here)
           end
           
       elseif PARAM.volume==0
           
           if PARAM.legendre==1
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [~,res(i,l)] = dropLegendreCurvilinearModes(t,xyMode,PARAM);
               cd(here)
           elseif PARAM.legendre==0
               here = pwd;
               cd('~/Documents/MATLAB/droplet_simulations/dropSpectral/')
               [~,res(i,l)] = dropExtensChebfunCurvilinearModes(t,xyMode,PARAM);
               cd(here)
           end
       end
       
   end
   
   %stop when there are no more data
   if i~=ite
   if T(k+1)>Tbreak||(T(k+1)==0&&k>1)
       ForceFree = ForceFree(1:i);
       V = V(1:i,l);
       t = T(1:i,l);
       res = res(1:i,l);
       display('Break')
       break;
   end
   end
    
end

figure(3)
if l>1
hold on
end
errV = (V(1:ite,l)-V0)/V0;
semilogy(t,abs(errV))
xlabel('t')
ylabel('err_V')
grid on
hold off
title('Error On Volume')
legend(cellLegend,'Location','Best')

figure(4)
if l>1
hold on
end
semilogy(t,abs(ForceFree(1:ite,l)))
xlabel('t')
ylabel('F')
title('Force On Drop')
grid on
hold off
legend(cellLegend,'Location','Best')

figure(5)
if l>1
hold on
end
semilogy(t,res(1:ite,l))
xlabel('t')
ylabel('res')
title('residuals')
grid on
hold off
legend(cellLegend,'Location','Best')

end













