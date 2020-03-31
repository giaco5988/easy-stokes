%validation extensional flow with stone

clear variables
close all

%paramters for upload
Ca1 = 0.01:0.01:0.16;
Ca2 = 0.01:0.01:0.12;
Ca3 = 0.01:0.01:0.1;
Ca = {Ca1 Ca2 Ca3};
lambda = [0.1 1 10];
dt = [5e-3 1e-2 5e-2];
manyTend = [100 100 400];
n = 50;
nameValX = {'xxx_01.mat' 'xxx_1.mat' 'xxx_10.mat'};
nameValY = {'yyy_01.mat' 'yyy_1.mat' 'yyy_10.mat'};

%directory
dir = '~/Documents/MATLAB/droplet_simulations/results/oneDropSpectral/validationExtensional/';
dirValidation = '~/Documents/MATLAB/droplet_simulations/results/extensional/viscosity';

for i = 1:numel(lambda)
    
    CaHere = Ca{i};
    D = zeros(numel(CaHere),1);
    
    for k = 1:numel(CaHere)
    
        %load data
        name = ['DropSpectral_ODE=' num2str(2) '_Legendre=' num2str(1) '_BC=' num2str(1) '_Ca=' num2str(CaHere(k)) '_visc=' num2str(lambda(i)) '_n=' num2str(n) '_D=' num2str(0) '_maxDT=' num2str(dt(i)) '_volCorr=' num2str(0) '_Tend=' num2str(manyTend(i)) '.mat'];
        load([dir name])
        
        %modes of last shape
        xMode = Y(end,1:2:end-1);
        yMode = Y(end,2:2:end);
        
        %from modes to points
        [x,y] = fromModesToGrid(xMode,yMode,PARAM);
        
        %compute deformation parameter
        L = max(x)-min(x);
        B = 2*max(y);
        D(k) = (L-B)/(L+B);
        
    end
    
    figure(1)
    if i>1
        hold on
    end
    plot(CaHere,D,'o-')
    xlabel('Ca')
    ylabel('D')
    title('Spectral Code Validation')
    grid on
    
end

for i = 1:numel(lambda)

    %upload results stone
    here = pwd;
    cd(dirValidation)
    upload = load(nameValX{i});
    xxx = upload.xxx;
    upload = load(nameValY{i});
    yyy = upload.yyy;
    cd(here)
    
    %plot results from stone
    hold on
    plot(xxx,yyy,'k.','MarkerSize',20)
    grid on
    
end

legend('\lambda=0.1','\lambda=1','\lambda=10','Stone et al 1989','Location','Best')











