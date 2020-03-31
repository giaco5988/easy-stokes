%check convergence on number of nodes

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%parameters
visc = 0.5;
Ca = 6;
amp = -1;
DELTA = 0.005;

%test on nodes
ELEM = [200 500];
LOOP = [10000 10000];

folder = ['~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=' num2str(visc) '_Ca=' num2str(Ca)];
here = pwd;

figure
hold on
for i = 1:numel(ELEM)
    
        cd(folder)

        load(['G=' num2str(amp) '_q=' num2str(ELEM(i)) '_visc=' num2str(visc) '_Dt=0.001_loop=' num2str(LOOP(i)) '_DELTA='...
            num2str(DELTA) '_Ca=' num2str(Ca) '_RK2.mat'])
        
        cd(here)
        
        %compute real area of the unperturbed sphere
        a = risa(:,1)';
        b = risb(:,1)';
        Volume = axis_int_gauss(a,b);
        R = nthroot(3/4/pi*Volume,3);
        A0 = 4*pi*R^2;
        
        %compute area variation norm
        all_area = Area-A0;
        Anorm = sqrt(all_area/all_area(1));
        
        %plot
        plot(0:checkpoint*deltaT:deltaT*loop,Anorm)
        xlabel('t')
        ylabel('||A||')
        grid on
        axis([0 6 1 1.5])
        title(['\delta=' num2str(DELTA)])
        
end

legend('q=200','q=500')