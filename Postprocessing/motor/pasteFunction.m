%paste simulations and get data

function [Displ,Vavg,Texit,DisplDIM,VavgDIM,TexitDIM] = pasteFunction(Ca,parts,path,theta,TimeEnd)

    set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

    %path = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/micromotor/parametric_study/theta005/';

    %simulations parameters
    %theta = -0.05;

    %which time to plot
    time = linspace(0,TimeEnd*parts,12);

    %load data for part 1
    if parts>1
        filename = [path 'Part1_LabFrame_mass=1_theta=' num2str(theta) '_el=320_dt=5e-05_visc=1_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=4_RK2.mat'];
    else
        filename = [path 'LabFrame_mass=1_theta=' num2str(theta) '_el=320_dt=0.0001_visc=1_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=4_RK2.mat'];
    end
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
    
    %dimensional
    radius = PARAM.R; %micromotor radius
    DisplDIM = radius*Displ;
    U = PARAM.mass/4/pi/radius^2;
    VavgDIM = Vavg*U;
    gamma = abs(U)*PARAM.visc2/PARAM.Ca;
    
    ind = find(RISa(PARAM.m+2,:)<0.5+RISa(PARAM.m+1,:),1,'first');
    %in case ind is empty take last instant
    if isempty(ind)==1
        ind = length(RISa);
        %ind = ;
    end
    Texit = PARAM.deltaT*PARAM.checkpoint*ind;
    
    %dimensional
    T = radius/gamma;
    TexitDIM = Texit*T;

    figure
    plot(ttt,v,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc)])

    hold on
    plot(ttt,Vavg*ones(1,ite),'--r','LineWidth',2)
    grid on
    xlabel('time')
    %axis([0 time(end) min(risy(m+q+2,:)) max(risy(m+q+2,:))])
    ylabel('v')
    legend('vel instant','vel average')

end