%plot DNS at some desired times

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

close all
clear variables

path = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/micromotor/parametric_study/';

%simulations parameters
Ca = 0.003;
theta = -0.05;

%plot snapshots
SubPlot = 1;

%when I have run many separate simulations
parts = 2;

%which time to plot
time = linspace(0,20*parts,12);

%figure in pixel
width = 1400;
height = 900;

%load data for part 1
filename = [path 'Part1_LabFrame_mass=1_theta=' num2str(theta) '_el=320_dt=5e-05_visc=1_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=4_RK2.mat'];
load(filename)
%velocity data
ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1)-1;
v = zeros(1,ite);

SIZE = size(risa);
SIZEy = size(risy);

RISa = zeros(SIZE(1)*2,SIZE(2)*parts-parts);
RISb = zeros(SIZE(1)*2,SIZE(2)*parts-parts);
RISy = zeros(SIZEy(1)*2,SIZEy(2)*parts-parts);

%paste different simulations
for k = 1:parts
    
    if k==1
        
        RISa(1:SIZE(1),1:ite/parts) = risa(:,1:end-1);
        RISb(1:SIZE(1),1:ite/parts) = risb(:,1:end-1);
        RISy(1:SIZEy(1),1:ite/parts) = risy;
        
    else
        
        %load data
        filename = [path 'Part' num2str(k) '_LabFrame_mass=1_theta=' num2str(theta) '_el=320_dt=5e-05_visc=1_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=4_RK2.mat'];
        load(filename)
        
        SIZE = size(risa);
        SIZEy = size(risy);
        
        RISa(1:SIZE(1),(1:ite/parts) + ite*(k-1)/parts) = risa(:,1:end-1);
        RISb(1:SIZE(1),(1:ite/parts) + ite*(k-1)/parts) = risb(:,1:end-1);
        RISy(1:SIZEy(1),(1:ite/parts) + ite*(k-1)/parts) = risy;
        
    end
    
end

for i = 1:ite
    %display(i)
    indVel = find(RISy(:,i)==0,2,'first');
    v(i) = RISy(indVel(2)-1,i);
end
ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end)-PARAM.checkpoint*PARAM.deltaT;

%compute velocity average
%integration
INT = ([diff(ttt),0]+[0,diff(ttt)])/2;
Displ = INT*v';
Vavg = 1/(ttt(end)-ttt(1))*INT*v';

figure(2)
plot(ttt,v,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc)])

hold on
plot(ttt,Vavg*ones(1,ite),'--r','LineWidth',2)
grid on
xlabel('time')
%axis([0 time(end) min(risy(m+q+2,:)) max(risy(m+q+2,:))])
ylabel('v')
legend('vel instant','vel average')

if SubPlot==1
for i = 1:numel(time)
    
    display(i)
    
    ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1);
    if i==numel(time)
        ite = ite-1;
    end
    
    %often used
    m = PARAM.m;    %q = nbr_el(ite);
    indNode = find(RISb(:,ite)==0,2,'first');
    q = indNode(2)-m-3;
    
    aMotor = RISa(1:m+1,ite);       bMotor = RISb(1:m+1,ite);
    aDrop = RISa(m+2:m+q+2,ite);    bDrop = RISb(m+2:m+q+2,ite);
    
    %velocity data
    vnow = v(1:ite);
    ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(i);
    
    figure(1)
    fig = figure(1);
    
    fig.Position = [400 200 width height];
    subplot(numel(time)/4,4,i)
    plot(aDrop,bDrop,'r',aDrop,-bDrop,'r')
    hold on
    plot(aMotor,bMotor,'b',aMotor,-bMotor,'b')
    hold off
    grid on
    axis equal
    axis([-4 12 -3 3])
    xlabel('x')
    ylabel('r')
    title(['t=' num2str(time(i))])
    
        
    hold on
    figure(2)
    plot(ttt(end),vnow(end),'.','MarkerSize',40)
    
end

%legend('v',['t=' num2str(time(1))],['t=' num2str(time(2))],['t=' num2str(time(3))],['t=' num2str(time(4))],'Location','Best')
title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta) ' x_0=' num2str(PARAM.start)])

end