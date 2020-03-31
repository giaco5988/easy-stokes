%study on the convergence of newton method

close all
clear variables

%path load data
%path = '~/Documents/MATLAB/droplet_simulations/results/rising_droplet/newton_method/convergence/';
path = '~/Documents/MATLAB/droplet_simulations/results/';

%physics data
Ca = 6;
lambda = 1;
elem = [20 30 40 50 60];
cellLegend = cell(size(elem));

figure
for i = 1:numel(elem)
    
    %store element in cell for legend
    cellLegend(i) = {num2str(elem(i))};

    %load data
    filename = [path 'RisingSpectralNewton_n=' num2str(elem(i)) '_Ca=' num2str(Ca) '_visc=' num2str(lambda) '.mat'];
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
title('Convergence rising droplet')
legend(cellLegend)