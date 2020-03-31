%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variables
Ca = '1.000000';
lambda = '1.000000';

%options
plotShape = 1;
plotCurv = 0;

%ID
D = '0.200000';
elem = '200';
totLoop = '500';
RK = '2';
CPUs = [1 2 4 8 12 16];

ManyArea = zeros(numel(CPUs),str2double(totLoop));
ManyV = zeros(numel(CPUs),str2double(totLoop));
SimulationTime = zeros(1,numel(CPUs));

%path
here = pwd;
path = '~/Documents/C++/rising_droplet/results/';

for i = 1:numel(CPUs)
    
    dt = '0.010000';
    
    results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' D '_elem=' elem '_dt=' dt '_loop=' totLoop '_RK=' RK '_CPUs=' num2str(CPUs(i)) ''];
    loop = str2double(totLoop);
    step = 1;

    display([num2str(i) ' of ' num2str(numel(CPUs))])
    
    cd(results)
    dt = str2double(dt);
    
    for k = 1:loop
        
        %display([num2str(k) ' of ' num2str(loop)])

        name = ['drop' num2str(k-1) '.txt'];

        A = importdata(name);

        x = A.data(:,1); y = A.data(:,2);

        %compute center of mass
        xcm = center_mass(x,y);

        %compute and save area
        ManyArea(i,k) = surf_gauss_vect(x',y');

        %compute and save volume
        ManyV(i,k) = axis_int_gauss_vect(x',y');
        
    end
    
    %get time
    name = ['time' num2str(k-1) '.txt'];
    A = importdata(name);
    SimulationTime(i) = A.data(1,1);
    
    %plot area variation norm and volume error
    R = nthroot(ManyV(i,1)/4/pi*3,3);
    A0 = 4*pi*R^2;
    myNorm = sqrt((ManyArea(i,:)-A0)/(ManyArea(i,1)-A0));
    figure(1)
    subplot(2,1,1)
    plot(0:dt:dt*(numel(myNorm)-1),myNorm)
    grid on
    ylabel('||A||')
    xlabel('t')
    hold on

    errV = (ManyV(i,:)-ManyV(1,1))/ManyV(i,1);
    subplot(2,1,2)
    plot(0:dt:dt*(numel(myNorm)-1),errV)
    grid on
    ylabel('err_V')
    xlabel('t')
    hold on
   
end

%compute speed up from simualtion time
SpeedUp = SimulationTime(1)./SimulationTime;

figure
plot(CPUs,SpeedUp,'o-')
grid on
hold on
plot([1 CPUs(end)],[1 CPUs(end)])
xlabel('CPUs')
ylabel('Speed Up')
legend('real','ideal','Location','Best')

cd(here)









