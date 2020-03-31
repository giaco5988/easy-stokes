%post processing of test time stepping spectral

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);
bemDir = '~/Documents/MATLAB/droplet_simulations/dropSpectral';
varDir = '~/Documents/MATLAB/droplet_simulations/results/edgeTracking/forPaper/variables';

%parameters
%manyCaUp = [.195 .175 .135 .105 .105]/2;    %n = 50
%manyVisc = [0.05 0.1 0.5 5 10];    %n = 50
manyCaUp = [0.225 .195 .175 .135 .105 .105]/2;    %n = 100
manyVisc = [0.02 0.05 0.1 0.5 5 10];    %n = 100
n = 100;
ODE = 2;        % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 1;         % 1 is extensional flow, 2 is rising droplet
volCorr = 0;    
manyMaxDT = [1e-3 1e-3 1e-3 2e-3 5e-3 1e-2];   %n = 100
%manyMaxDT = [2e-3 2e-3 2e-3 2e-2 2e-2];    %n = 50
legendre = 1;
Dup = 0;          % initial shape
CPUs = 16;
res = 2;        %results or server
convergeShape = 0.15;

%n=100
TendCell{1} = 255:265;
TendCell{2} = 215:225;
TendCell{3} = 290:310;
TendCell{4} = 290:305;
TendCell{5} = 1150:1165;
TendCell{6} = 1270:1275;

%n=50
% TendCell{1} = 115:125;
% TendCell{2} = 290:310;
% TendCell{3} = 295:310;
% TendCell{4} = 1520:5:1570;
% TendCell{5} = 1250:5:1290;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==2
    dir = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/edgeTracking/forPaper/dataElongationExperiment/';
end

%initialize
Lavg = zeros(numel(manyVisc),1);
LerrUp = zeros(numel(manyVisc),1);
LerrDown = zeros(numel(manyVisc),1);

Lbefore = 0;
for i = 1:numel(manyCaUp)
    
    display([num2str(i) ' of ' num2str(numel(manyCaUp))])
    
    Ca = manyCaUp(i);
    visc = manyVisc(i);
    Tplot = TendCell{i};
    maxDT = manyMaxDT(i);
        
    first = 1;
        for l = 1:numel(Tplot)
            
            TendUp = Tplot(l);
        
            %filename
            name = ['DropSpectral_ODE=' num2str(ODE) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '_Ca=' num2str(Ca) '_visc=' num2str(visc) '_n=' num2str(n) '_D=' num2str(Dup) '_maxDT=' num2str(maxDT) '_volCorr=' num2str(volCorr) '_Tend=' num2str(TendUp) '.mat'];

            %upload data
            upload = [dir name];
            load(upload)

            %take first shape
            xMode = Y(1,1:2:end-1);
            yMode = Y(1,2:2:end);
            x = fromModesToGrid(xMode,yMode,PARAM);

            %elongation of he initial shape
            Lstone = max(x);

            %take last shape
            %ind = find(t==0,2,'first');
            xMode = Y(end,1:2:end-1);
            yMode = Y(end,2:2:end);
            [x,y] = fromModesToGrid(xMode,yMode,PARAM);

            %compute deformation parameter of last shape
            L = 2*max(x);
            B = 2*y(round(numel(y)/2));
            D = (L-B)/(L+B);

            if PARAM.BC==1
                here = pwd;
                cd(bemDir)
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
                cd(here)
                [maxCa,ind] = max(manyCa);
                manyD = manyD(1:ind);
                manyCa = manyCa(1:ind);
                [~,ind] = min(abs(PARAM.Ca-manyCa));
                Dca = manyD(ind);
                OutIn  = abs(D-Dca)/Dca > convergeShape;

            elseif PARAM.BC==2

                OutIn  = abs(D) > convergeShape;

            end

            if OutIn==1
               symbol = 'xr'; 
            elseif OutIn==0
               symbol = 'ok';
            end

            semilogx(visc,Lstone,symbol)
            hold on
            xlabel('\lambda')
            ylabel('L')
            grid on
            
            if OutIn==1 && first==1
               Lavg(i) = (Lstone+Lbefore)/2;
               LerrUp(i) = Lstone;
               LerrDown(i) = Lbefore;
               first = 0;
            end
            Lbefore = Lstone;
        
        end
    
end

here = pwd;
cd(varDir)
load('saddleNode.mat')
load('manyLambdaBif.mat')
load('Lup.mat')
load('Ldown.mat')
cd(here)

semilogx(manyLambda,Lup)

figure
semilogx(manyVisc,Lavg)
hold on
semilogx(manyLambda,Lup)
semilogx(manyLambda,Ldown)
%semilogx(manyLambda,LsaddleNode)
%semilogx(manyVisc,LerrUp,'kx')
%semilogx(manyVisc,LerrDown,'kx')
%plot error bar
for i =1 :numel(manyVisc)
    semilogx(ones(2,1)*manyVisc(i),[LerrDown(i) LerrUp(i)],'o-k');
end
xlabel('\lambda')
legend('DNS, step reduction of Ca','Edge state elongation','Stable shape elongation','Location','Best')
ylabel('L')
grid on
title('Droplet Elongation')

figure
if n==100
    errLup = abs(Lup([1 2 3 4 6 7])-Lavg)./Lavg;
    errLdown = abs(Ldown([1 2 3 4 6 7])-Lavg)./Lavg;
elseif n==50
    errLup = abs(Lup([2 3 4 6 7])-Lavg)./Lavg;
    errLdown = abs(Ldown([2 3 4 6 7])-Lavg)./Lavg;
end

semilogx(manyVisc,errLup)
%loglog(manyVisc,LerrUp,'x')
hold on
semilogx(manyVisc,errLdown)
%loglog(manyVisc,LerrDown,'x')
%loglog(manyVisc,Lavg)
%loglog(manyLambda,Lup)
xlabel('\lambda')
legend('Error compared to the edge state','Error compared ro the stable staedy state','Location','Best')
%axis([1e-2 10 1e-1 1])
ylabel('\Delta L/L')
grid on






