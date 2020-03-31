%plot DNS at some desired times

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

close all
clear variables

theta = -0.05;
TimeEnd = 24;
path = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/micromotor/parametric_study/theta005/';

%simulations parameters
%Ca = 0.001:0.001:0.006;
Ca = 0.001:0.001:0.006;

%initialize
Displ = zeros(1,numel(Ca));
Vavg = zeros(1,numel(Ca));
Texit = zeros(1,numel(Ca));
DisplDIM = zeros(1,numel(Ca));
VavgDIM = zeros(1,numel(Ca));
TexitDIM = zeros(1,numel(Ca));

for k = 1:numel(Ca)
    
    display([num2str(k) ' of ' num2str(numel(Ca))])
   
    if theta==-0.05
        
%         if Ca(k)>0.003
%             parts = 3;
%         else
%             parts = 2;
%         end
        parts = 1;
        
    elseif theta==-0.02
        
       parts = 1; 
        
    end
    
    [Displ(k),Vavg(k),Texit(k),DisplDIM(k),VavgDIM(k),TexitDIM(k)] = pasteFunction(Ca(k),parts,path,theta,TimeEnd);
    
end

figure
plot(Ca,Displ,'ok-')
xlabel('Ca')
ylabel('displace')
grid on

figure
plot(Ca,Vavg,'ok-')
xlabel('Ca')
ylabel('V_{avg}')
grid on

figure
plot(Ca,Texit,'o-k')
xlabel('Ca')
ylabel('T_{exit}')
grid on

figure
plot(Ca,DisplDIM,'or-')
xlabel('Ca')
ylabel('displace [m]')
grid on

figure
plot(Ca,VavgDIM,'or-')
xlabel('Ca')
ylabel('V_{avg} [m/s]')
grid on

figure
plot(Ca,TexitDIM,'o-r')
xlabel('Ca')
ylabel('T_{exit} [s]')
grid on