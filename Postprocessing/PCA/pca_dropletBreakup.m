clear variables
close all

%load data
num_modes = 10;
Kproject = 3;
load(['/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/data_' num2str(num_modes) 'Modes_CM_trajectories.mat'])

%number of examples
[m, ~] = size(data);

%features and labels
X = data(:, 1:end-1);
y = data(:, end);

%normalize
%X = (X-mean(X))./std(X);

%covariance matrix
cov = 1/m*(X'*X);

%compute eigenvalues
[V,D] = eig(cov);

%project data
Z = X*V(:,1:Kproject);
    
figure
plot(X(y==1,1),X(y==1,2),'.')
grid on
xyLabelTex('f_1','f_2')
title('Phase space, first 2 Legendre')
hold on
plot(X(y==0,1),X(y==0,2),'.')
legend('Breakup','No breakup','Location','Best')
mu = mean(X);
drawLine(mu, mu + 5*D(1,1) * V(:,1)', '-k', 'LineWidth', 2);
drawLine(mu, mu + 5*D(2,2) * V(:,2)', '-k', 'LineWidth', 2);
    
if Kproject==2
    
    figure
    plot(Z(y==1,1),Z(y==1,2),'.')
    hold on
    plot(Z(y==0,1),Z(y==0,2),'.')
    legend('Breakup','No breakup','Location','Best')
    grid on
    xyLabelTex('V_1','V_2')
    title('PCA projection')

elseif Kproject==3
    
    figure
    plot3(Z(y==1,1),Z(y==1,2),Z(y==1,3),'.')
    hold on
    plot3(Z(y==0,1),Z(y==0,2),Z(y==0,3),'.')
    legend('Breakup','No breakup','Location','Best')
    grid on
    xyzLabelTex('V_1','V_2','V_3')
    title('PCA projection')
    
end
   
   
   
   
   