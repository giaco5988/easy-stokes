%plot droplet relaxation

close all
clear all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

load('~/Documents/MATLAB/droplet_simulations/results/relaxation/Relaxation_q=200_visc=1_Dt=0.1_loop=3000_DELTA=5_Ca=10_RK2.mat')

%time at which plotting
time = [0:25:100 200];

for i = 1:numel(time)
    
    disp(i)
    
    ite = time(i)/deltaT;
    if ite==0
        ite=1;
    end
    a = risa(:,ite);    b = risb(:,ite);
    
    shape_x = [a; -a];   shape_y = [b; -b];
    
    hold all
    plot(shape_x,shape_y)
    axis equal
    %hold off
    
end

grid on
xlabel('x')
ylabel('r')
legend('t=0','t=2.5','t=5','t=7.5','t=10','t=20','Location','Best')
axis([-3.5 3.5 -1.5 1.5])