%uplaod data drom drop from c++

clear variables
close all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%physical variables
Ca = '0.050000';
lambda = '1.000000';

%options
plotShape = 0;
plotCurv = 0;

%ID
D = '0.020000';
MANYelem = [500 1000 2000 3000 4000 5000];
totLoop = '1';
RK = '1';
CPUs = [1 2];
checkpoint = 1;

figure(1)
plot([1 CPUs(end)],[1 CPUs(end)])
hold on

figure(2)
plot([1 CPUs(end)],[1 CPUs(end)])
hold on

ManyArea = zeros(numel(CPUs),str2double(totLoop)/checkpoint);
ManyV = zeros(numel(CPUs),str2double(totLoop)/checkpoint);
MatrixTime = zeros(1,numel(CPUs));
SolveTime = zeros(1,numel(CPUs));

%path
here = pwd;
path = '~/Documents/C++/drop_extensionalMPI/results/';

for l = 1:numel(MANYelem)
    
    elem = MANYelem(l);
    elem = num2str(elem);

    for i = 1:numel(CPUs)

        dt = '0.010000';

        results = [path 'Ca=' Ca '_lambda=' lambda '_delta=' D '_elem=' elem '_dt=' dt '_loop=' totLoop '_RK=' RK '_CPUs=' num2str(CPUs(i)) ''];
        loop = str2double(totLoop)/checkpoint;
        step = 1;

        display([num2str(i) ' of ' num2str(numel(CPUs))])

        cd(results)
        dt = str2double(dt);

        %get time matrix building
        name = 'BuildMatrix.txt';
        A = importdata(name);
        MatrixTime(i) = A(1,1);

        %get time matrix building
        name = 'SolveMatrix.txt';
        A = importdata(name);
        SolveTime(i) = A(1,1);

    end

%compute speed up from simualtion time
SpeedUpMatrix = MatrixTime(1)./MatrixTime;
SpeedUpSolve = SolveTime(1)./SolveTime;

figure(1)
grid on
plot(CPUs,SpeedUpMatrix,'o-')
xlabel('CPUs')
ylabel('Speed Up')
legend('ideal','500','1000','2000','3000','4000','5000','Location','northwest')
title('SpeedUP matrix building')

figure(2)
grid on
plot(CPUs,SpeedUpSolve,'o-')
xlabel('CPUs')
ylabel('Speed Up')
legend('ideal','500','1000','2000','3000','4000','5000','Location','northwest')
title('SpeedUP matrix solving')

figure(3)
grid on
plot(CPUs,MatrixTime,'o-')
hold on
xlabel('CPUs')
ylabel('Time [ms]')
legend('500','1000','2000','3000','4000','5000','Location','northeast')
title('Time matrix building')

figure(4)
grid on
plot(CPUs,SolveTime,'o-')
hold on
xlabel('CPUs')
ylabel('Time [ms]')
legend('500','1000','2000','3000','4000','5000','Location','northeast')
title('Time matrix solving')

cd(here)

end









