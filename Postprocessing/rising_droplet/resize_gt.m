%resizing G(t) HAS TO BE USED IN ~/Documents/MATLAB/droplet_simulations/results/delta_graph/data

clear all
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%delta that I want to analyze (before the stables and after the unstables)
delta = [-0.0036 0.00375 -0.004 0.0058 0.0063 -0.0077 -0.0037 0.004 -0.00425 0.0059 0.0065 -0.0078];
%time at wchich the norm maximizes the gain
time = [0.25 0.55 1.15 3.3 5.35 6.1 0.25 0.55 1.15 3.3 5.35 6.1];
%when in the curve: before, maximum or after
when = [1 1 1 2 3 3];
GT = [0.1 0.25 0.5 1 0.5 0.1];

special = [0.002874 0.002948 0.002574 0.002106 0.001541 0.002141];
A_trunc0 = 4*pi+2*pi*special;

%preallocation
V_stab_lin = zeros(numel(delta)/2,1);
A_stab_lin = zeros(numel(delta)/2,1);
A0_stab_lin = zeros(numel(delta)/2,1);
V_stab_nonlin = zeros(numel(delta)/2,1);
V_unstab_lin = zeros(numel(delta)/2,1);
A_unstab_lin = zeros(numel(delta)/2,1);
A0_unstab_lin = zeros(numel(delta)/2,1);
V_unstab_nonlin = zeros(numel(delta)/2,1);
A_trunc = zeros(numel(delta)/2,1);
%A0_trunc = zeros(numel(delta)/2,1);

for i = 1:numel(delta)/2
    
    disp(i)
    
    filename = strcat('D=',num2str(delta(i)));
    cd(filename)
    
    filename = strcat('crd_rzlam0.5_ca6_mA1000_delta',num2str(delta(i)),'_time',num2str(time(i)),'.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    V1 = axis_int_gauss(y',x');
    A_stab_lin(i) = surf_gauss(y',x');
    
    filename = strcat('crd_rzlam0.5_ca6_mA1000_delta',num2str(delta(i)),'_time0.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    V2 = axis_int_gauss(y',x');
    A0_stab_lin(i) = surf_gauss(y',x');
    
    V_stab_lin(i) = V1/V2;
    
    cd ..
    filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
    cd(filename)
    
    filename = strcat('crd_rzlam0.5_ca6_mA1000_delta',num2str(delta(i+numel(delta)/2)),'_time',num2str(time(i)),'.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    V1 = axis_int_gauss(y',x');
    A_unstab_lin(i) = surf_gauss(y',x');
    
    filename = strcat('gt_tlam0.5_ca6_mA1000_delta',num2str(delta(i+numel(delta)/2)),'.dat');
    A = importdata(filename);
    
    t = A(:,1);
    area_norm = A(:,2);
    
    [~,l] = min(abs(time(i)-t));
    
    A_trunc(i) = 4*pi+area_norm(l).^2*2*pi*special(i);
    
    filename = strcat('crd_rzlam0.5_ca6_mA1000_delta',num2str(delta(i+numel(delta)/2)),'_time0.dat');
    A = importdata(filename);
    
    x = A(:,1);
    y = A(:,2);
    
    V2 = axis_int_gauss(y',x');
    A0_unstab_lin(i) = surf_gauss(y',x');
    
    V_unstab_lin(i) = V1/V2;
    
    if when(i)==1
        
        if i==1
        
            cd ..
            filename = strcat('D=',num2str(delta(i)));
            cd(filename)

            t = [0 0.0625 0.125 0.1875 0.25 1.22 2.19 3.16 4.13 5.1 6.07 7.04 8.01 8.98 9.95];
            plot_area_norm(t,delta(i));
            
            cd ..
            filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
            cd(filename)
            plot_area_norm(t,delta(i+numel(delta)/2));
            
        elseif i==2
            
            cd ..
            filename = strcat('D=',num2str(delta(i)));
            cd(filename)

            t = [0 0.1375 0.275 0.55 1.49 2.43 3.37 4.31 5.25 6.19 7.13 8.07 9.01 9.95];
            plot_area_norm(t,delta(i));
            
            cd ..
            filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
            cd(filename)
            plot_area_norm(t,delta(i+numel(delta)/2));
            
        elseif i==3
            
            cd ..
            filename = strcat('D=',num2str(delta(i)));
            cd(filename)

            t = [0 0.2875 0.575 0.8625 1.15 2.03 2.91 3.79 4.67 5.55 6.43 7.31 8.19 9.07 9.95];
            plot_area_norm(t,delta(i));
            
            cd ..
            filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
            cd(filename)
            plot_area_norm(t,delta(i+numel(delta)/2));
            
        end
        
        cd ..
        cd ..
        cd before
        filename = strcat('before_Gt=',num2str(GT(i)),'_q=200_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i)),'_Ca=6_RK2.mat');
        load(filename)
        
        V_stab_nonlin(i) = V(round(time(i)/deltaT));
        
        filename = strcat('before_Gt=',num2str(GT(i)),'_q=0_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i+numel(delta)/2)),'_Ca=6_RK2.mat');
        load(filename)
        cd ..
        cd data
        
        V_unstab_nonlin(i) = V(round(time(i)/deltaT));
        
    elseif when(i)==2
        
        cd ..
        filename = strcat('D=',num2str(delta(i)));
        cd(filename)
        
        t = [0 0.825 1.65 2.475 3.3 4.125 4.95 5.775 6.6 7.425 8.25 9.075 9.9 10.725];
        plot_area_norm(t,delta(i));
        
        cd ..
        filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
        cd(filename)
        plot_area_norm(t,delta(i+numel(delta)/2));
        
        cd ..
        cd ..
        cd max
        filename = strcat('max_q=200_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i)),'_Ca=6_RK2.mat');
        load(filename)
        
        V_stab_nonlin(i) = V(round(time(i)/deltaT));
        
        filename = strcat('max_q=0_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i+numel(delta)/2)),'_Ca=6_RK2.mat');
        load(filename)
        cd ..
        cd data
        
        V_unstab_nonlin(i) = V(round(time(i)/deltaT));
        
     elseif when(i)==3
         
        cd ..
        filename = strcat('D=',num2str(delta(i)));
        cd(filename)
        
        if i==5
            
            t = [0 1.3375 2.675 4.0125 5.35 6.6875 8.025 9.3625 10.7];
            plot_area_norm(t,delta(i));
            
            cd ..
            filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
            cd(filename)
            plot_area_norm(t,delta(i+numel(delta)/2));
            
        elseif i==6
        
            t = [0 1.525 3.05 4.575 6.1 7.625 9.15 10.675];
            plot_area_norm(t,delta(i));
            
            cd ..
            filename = strcat('D=',num2str(delta(i+numel(delta)/2)));
            cd(filename)
            plot_area_norm(t,delta(i+numel(delta)/2));
            
        end
        
        cd ..
        cd ..
        cd after
        filename = strcat('after_Gt=',num2str(GT(i)),'_q=200_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i)),'_Ca=6_RK2.mat');
        load(filename)
        
        V_stab_nonlin(i) = V(round(time(i)/deltaT));
        
        filename = strcat('after_Gt=',num2str(GT(i)),'_q=0_visc=0.5_Dt=0.01_loop=1500_DELTA=',num2str(delta(i+numel(delta)/2)),'_Ca=6_RK2.mat');
        load(filename)
        cd ..
        cd data
        
        V_unstab_nonlin(i) = V(round(time(i)/deltaT));
        
    end
            
end

%plot volumes for linear simulations
% figure
% plot(time(1:numel(time)/2),V_stab_lin,time(1:numel(time)/2),V_unstab_lin)
% xlabel('time')
% ylabel('Volume in linear simulations')
% legend('Less stable','Less unstable')

%plot G(t) as if V_nl=1
filename = 'Gt_area_t_Ca6_lam0.5_mA1000.dat';
A = importdata(filename);
t = A(:,1);
g = A(:,2);
%gt = zeros(numel(delta)/2,1);

% for i=1:numel(time)/2
%     
%     [~,k] = min(abs(time(i)-t));
%     gt(i) = g(k);
%     
% end

A_trunc_res = A_trunc.*(V_unstab_nonlin./V_unstab_lin).^(2/3);
A_norm_trunc_res = sqrt((A_trunc_res-4*pi)'./(A_trunc0-4*pi));

% I MISS SPECIAL FOR THE TIME0!!!
%gt_special_res = ;

%plot G(t) fr0m linear simulations data with complete area
figure
plot(t,g,time(1:numel(time)/2),sqrt((A_stab_lin-4*pi)./(A0_stab_lin-4*pi)),'o-',time(1:numel(time)/2),A_norm_trunc_res,'o-m')
xlabel('time')
ylabel('A norm')
grid on

resize = A_stab_lin.*(V_stab_nonlin./V_stab_lin).^(2/3);
gt_stab_res = sqrt((resize-4*pi)./(A0_stab_lin-4*pi));

resize = A_unstab_lin.*(V_unstab_nonlin./V_unstab_lin).^(2/3);
gt_unstab_res = sqrt((resize-4*pi)./(A0_unstab_lin-4*pi));

figure
%plot(t,g,time(1:numel(time)/2),gt,time(1:numel(time)/2),gt_stab_res,time(1:numel(time)/2),gt_unstab_res)
plot(t,g,time(1:numel(time)/2),gt_stab_res,'r',time(1:numel(time)/2),gt_unstab_res,'--')
xlabel('time')
ylabel('Area-variation norm')
legend('G(t)','G(t) stable resized','G(t) unstable resized','Location','Best')
grid on