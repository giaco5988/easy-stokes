%study on the convergence of newton method

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%path load data
%path = '~/Documents/MATLAB/droplet_simulations/results/extensional/newton_method/convergence/';
path = '~/Documents/MATLAB/droplet_simulations/results/';

%physics data
Ca = 0.1;
lambda = 1;
elem = [20 30 40 50 60];
cellLegend = cell(size(elem));

figure
for i = 1:numel(elem)
    
    %store element in cell for legend
    cellLegend(i) = {num2str(elem(i))};

    %load data
    filename = [path 'newtonMethodSpectralXYmodes_n=' num2str(elem(i)) '_Ca=' num2str(Ca) '_visc=' num2str(lambda) '.mat'];
    load(filename)
    
    %plot residuals
    semilogy(manyRES,'-o')
    if i==1
        hold on
    end
    
end
grid on
xlabel('iteration')
ylabel('R')
title('Convergence')
legend(cellLegend,'Location','Best')