%plot DNS at some desired times

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/';
source = '~/Documents/MATLAB/droplet_simulations/server/';

%choose what to compare
CompareVisc = 0;
CompareLenght = 0;
CompareCa = 0;
CompareDt = 0;
CompareElem = 0;
CompareHrep = 0;
CompareIntensRep = 1;
CompareResize = 0;

%real physical parameters (from Manjare et al. number 9 in TAB)
%RealGamma = [33*1e-3 33*1e-3];
RealGamma = 33*1e-3;
RealRadius = 2*1e-6;
RealVisc =8.9*1e-4;
InjectionRate  = 1.3*1e-14/1.2;

%compute experimental Ca
CaExperim = InjectionRate*RealVisc/RealRadius^2./RealGamma;
display(['Experimental Ca is equal to ' num2str(CaExperim)])

%parameters
if CompareVisc==1
    Ca = 0.01;
    inflate = 0;
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = 228;
    visc = [0.01 0.1 1];
    InPos = 4;
    dt = 0.01;
    Lenght = 10;
    time = 1550*ones(numel(visc),1);
    RK = 2;
elseif CompareLenght==1
    Ca = 0.001;
    inflate = 0;
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = [220 280 320];
    visc = 1;
    InPos = [2 3 4];
    dt = 0.0001;
    Lenght = [5 8 10];
    time = [10 13 16];
    RK = 2;
elseif CompareCa==1
    Ca = [0.001 0.01];
    inflate = 0;
    theta = -0.02;
    alpha = 0.8;
    elem = 228;
    visc = 0.1;
    InPos = 4;
    dt = 0.01;
    Lenght = 10;
    time = [5000 1300];
    RK = 2;
elseif CompareDt==1
    Resize = 1;
    MotorFree = 1;
    Ca = 0.01;
    inflate = 0;
    dt = [0.002 0.005 0.01];
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = 308;
    %elem = 178;
    visc = 0.1;
    Hrep = 3;
    IntensRep = 1e-6;
    InPos = 1;
    Lenght = 10;
    time = [390 900 900];
    RK = 2;
elseif CompareElem==1
    Resize = 1;
    MotorFree = 1;
    Ca = 0.01;
    inflate = 0;
    dt = 0.01;
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = [228 238 268 288 308 358];
    %elem = 178;
    visc = 0.1;
    Hrep = 3;
    IntensRep = 1e-6;
    InPos = 1;
    Lenght = 10;
    time = 1000*ones(numel(elem),1);
    RK = 2;
elseif CompareHrep==1
    Resize = 0;
    MotorFree = 1;
    Ca = 0.01;
    inflate = 0;
    dt = 0.01;
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = 188;
    visc = 1;
    InPos = 1;
    Lenght = 10;
    Hrep = [3 4];
    time = 1000*ones(numel(Hrep),1);
    RK = 2;
    IntensRep = 1e-4;
elseif CompareResize==1
    Resize = [0 1];
    MotorFree = 1;
    Ca = 0.01;
    inflate = 0;
    dt = 0.01;
    theta = -0.02;
    alpha = 0.8;
    mass = 1;
    elem = 188;
    visc = 1;
    InPos = 1;
    Lenght = 10;
    Hrep = 4;
    time = 1000*ones(numel(Hrep),1);
    RK = 2;
    IntensRep = 1e-4;
elseif CompareIntensRep==1
    deflationWall = 1;
    MotorFree = 0;
    Resize = 0;
    Ca = 0.01;
    inflate = 0;
    dt = 0.01;
    theta = -0.05;
    alpha = 0.8;
    mass = 1;
    elem = 308*ones(8,1);
    visc = 0;
    InPos = 1;
    Lenght = 10;
    Hrep = 6;
    repLenght = [0.1*ones(4,1); 0.01*ones(4,1)];
    IntensRep = -[100 10 1 0.1 100 10 1];
    time = 990*ones(numel(IntensRep),1);
    RK = 2;
    maxVel = 4e-3;      %axis for plotting
end

if CompareCa==1

    %loop in visc
    AVG = zeros(1,numel(Ca));
    AVGdim = zeros(1,numel(Ca));
    Vmax = zeros(1,numel(Ca));
    VmaxDim = zeros(1,numel(Ca));
    for k = 1:numel(Ca)
        
        %store element in cell for legend
        cellLegend(k) = {['Ca=' num2str(Ca(k))]};
        
        display([num2str(k) 'of ' num2str(numel(Ca))])

        load([source 'ConicalMotor_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc(1)) '_Ca=' num2str(Ca(k)) '_R=' num2str(1) '_L=' num2str(Lenght) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %reference scales to make it dimensional
        vScale = RealGamma(k)/RealVisc;
        %tScale = RealRadius*RealVisc/PARAM.Ca/RealGamma;
        tScale = RealRadius/vScale;
        %tScale = RealGamma/RealVisc;
        %vScale = RealRadius/tScale;
        
        %velocity data
        ite = round(time(k)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        for i = 1:ite
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(k);

        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Micromotor velocity \lambda=' num2str(visc)])
        hold on
        
        figure(2)
        plot(ttt*tScale,v*vScale)
        grid on
        xlabel('t[s]')
        ylabel('v[m/s]')
        title('Micromotor velocity')
        hold on
        
        %compute average velocity
        AVG(k) = (1/time(k))*trapz(ttt,v);
        AVGdim(k) = (1/time(k))*trapz(ttt,v)*vScale;
        
        %compute max velocity
        Vmax(k) = max(v);
        VmaxDim(k) = max(v)*vScale;

    end

    xi = Lenght/2/PARAM.R;
    
    figure(1)
    legend(cellLegend,'Location','Best')
    
    figure(2)
    legend(cellLegend,'Location','Best')
    
    %plot max velocity as a function of the capillary number
    figure
    plot(Ca,AVG,'o-')
    grid on
    xlabel('Ca')
    ylabel('V_{avg}')
    title('Average velocity')
    
    %plot max velocity as a function of the capillary number
    figure
    plot(Ca,AVGdim,'o-')
    grid on
    xlabel('Ca')
    ylabel('V_{avg}[m/s]')
    title('Average velocity')
    
    %plot max velocity as a function of the capillary number
    figure
    plot(Ca,Vmax,'o-')
    grid on
    xlabel('Ca')
    ylabel('V_{max}')
    title('Max velocity')
    
    %plot max velocity as a function of the capillary number
    figure
    plot(Ca,VmaxDim,'o-')
    grid on
    xlabel('Ca')
    ylabel('V_{max}[m/s]')
    title('Max velocity')

end

if CompareDt==1
    
    cellLegend = cell(size(CompareDt));

    %loop in visc
    AVG = zeros(1,numel(Ca));
    AVGdim = zeros(1,numel(Ca));
    Vmax = zeros(1,numel(Ca));
    VmaxDim = zeros(1,numel(Ca));
    for k = 1:numel(dt)
        
        %store element in cell for legend
        cellLegend(k) = {['dt=' num2str(dt(k))]};
        
        display([num2str(k) 'of ' num2str(numel(dt))])

        load([source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt(k)) '_visc=' num2str(visc(1)) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(Lenght) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %reference scales to make it dimensional
        vScale = PARAM.Ca*RealGamma/RealVisc;
        %tScale = RealRadius*RealVisc/PARAM.Ca/RealGamma;
        tScale = RealRadius/vScale;
        %tScale = RealGamma/RealVisc;
        %vScale = RealRadius/tScale;
        
        %velocity data
        ite = round(time(k)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        for i = 1:ite
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(k);

        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title('Micromotor velocity')
        hold on
        
        figure(2)
        plot(ttt*tScale,v*vScale)
        grid on
        xlabel('t[s]')
        ylabel('v[m/s]')
        title('Micromotor velocity')
        hold on
        
        %compute average velocity
        AVG(k) = (1/time(k))*trapz(ttt,v);
        AVGdim(k) = (1/time(k))*trapz(ttt,v)*vScale;
        
        %compute max velocity
        Vmax(k) = max(v);
        VmaxDim(k) = max(v)*vScale;

    end

    xi = Lenght/2/PARAM.R;
    
    figure(1)
    legend(cellLegend,'Location','Best')
    
    
    figure(2)
    legend(cellLegend,'Location','Best')

end

if CompareElem==1
    
    cellLegend = cell(size(CompareDt));

    %loop in visc
    AVG = zeros(1,numel(Ca));
    AVGdim = zeros(1,numel(Ca));
    Vmax = zeros(1,numel(Ca));
    VmaxDim = zeros(1,numel(Ca));
    for k = 1:numel(elem)
        
        display([num2str(k) 'of ' num2str(numel(elem))])

        load([source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem(k)) '_dt=' num2str(dt(1)) '_visc=' num2str(visc(1)) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(Lenght) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %store element in cell for legend
        cellLegend(k) = {['wall=' num2str(PARAM.m) ' drop=' num2str(PARAM.q)]};
        
        %reference scales to make it dimensional
        vScale = PARAM.Ca*RealGamma/RealVisc;
        %tScale = RealRadius*RealVisc/PARAM.Ca/RealGamma;
        tScale = RealRadius/vScale;
        %tScale = RealGamma/RealVisc;
        %vScale = RealRadius/tScale;
        
        %velocity data
        ite = round(time(k)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        ForceDrop = zeros(1,ite);
        for i = 1:ite
            
            %useful variables
            m = find(risa(2:end,i)==risa(1,i));
            indVel = find(risy(:,i)==0,2,'first');
            indXY = find(risb(:,i)==0,2,'first');
            aDrop = risa(m+2:indXY(2)-1,i);
            bDrop = risb(m+2:indXY(2)-1,i);
            
            %motor velocity
            v(i) = risy(indVel(2)-1,i);
            
            %check force free condition
            %ForceDrop(i) = forceOnDrop(aDrop,bDrop,KKK(1:nbr_el(i)+1,i),1);
        
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(k);

        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Micromotor velocity Ca=' num2str(Ca)])
        hold on
        
        figure(2)
        plot(ttt,ForceDrop)
        grid on
        xlabel('t')
        ylabel('F')
        title(['Force on drop Ca=' num2str(Ca)])
        hold on
        
%         figure(2)
%         plot(ttt*tScale,v*vScale)
%         grid on
%         xlabel('t[s]')
%         ylabel('v[m/s]')
%         title('Micromotor velocity')
%         hold on
        
        %compute average velocity
        AVG(k) = (1/time(k))*trapz(ttt,v);
        AVGdim(k) = (1/time(k))*trapz(ttt,v)*vScale;
        
        %compute max velocity
        Vmax(k) = max(v);
        VmaxDim(k) = max(v)*vScale;

    end

    xi = Lenght/2/PARAM.R;
    
    figure(1)
    legend(cellLegend,'Location','Best')
    
    figure(2)
    legend(cellLegend,'Location','Best')
    
    
%     figure(2)
%     legend(cellLegend,'Location','Best')

end

if CompareLenght==1

    %loop in visc
    AVG = zeros(1,numel(Lenght));
    AVGdim = zeros(1,numel(Lenght));
    Vmax = zeros(1,numel(Ca));
    VmaxDim = zeros(1,numel(Ca));
    for k = 1:numel(Lenght)
        
        display([num2str(k) 'of ' num2str(numel(Lenght))])

        load(['LabFrame_mass=' num2str(mass) '_theta=' num2str(theta) '_el=' num2str(elem(k)) '_dt=' num2str(dt) '_visc=' num2str(visc(1)) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(Lenght(k)) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos(k)) '_RK=' num2str(RK) '.mat']);

        %reference scales to make it dimensional
        vScale = PARAM.Ca*RealGamma/RealVisc;
        %tScale = RealRadius*RealVisc/PARAM.Ca/RealGamma;
        tScale = RealRadius/vScale;
        %tScale = RealGamma/RealVisc;
        %vScale = RealRadius/tScale;
        
        %velocity data
        ite = round(time(k)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        for i = 1:ite
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(k);

        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Micromotor velocity Ca=' num2str(Ca)])
        hold on
        
        figure(2)
        plot(ttt*tScale,v*vScale)
        grid on
        xlabel('t[s]')
        ylabel('v[m/s]')
        title(['Micromotor velocity Ca=' num2str(Ca)])
        hold on
        
        %compute average velocity
        AVG(k) = (1/time(k))*trapz(ttt,v);
        AVGdim(k) = (1/time(k))*trapz(ttt,v)*vScale;
        
        %compute max velocity
        Vmax(k) = max(v);
        VmaxDim(k) = max(v)*vScale;

    end

    xi = Lenght/2/PARAM.R;
    
    figure(1)
    legend(['\xi=' num2str(Lenght(1))],['\xi=' num2str(Lenght(2))],['\xi=' num2str(Lenght(3))],'Location','Best')
    %legend(['\xi=' num2str(xi(1))],['\xi=' num2str(xi(2))],'Location','Best')
    
    figure(2)
    legend(['\xi=' num2str(Lenght(1))],['\xi=' num2str(Lenght(2))],['\xi=' num2str(Lenght(3))],'Location','Best')
    
    %plot average velocity as a function of the geometry
    figure
    plot(xi,AVG,'o-')
    grid on
    xlabel('\xi')
    ylabel('V_{avg}')
    title('Average velocity')
    
    %plot average velocity as a function of the geometry
    figure
    plot(xi,AVGdim,'o-')
    grid on
    xlabel('\xi')
    ylabel('V_{avg}[m/s]')
    title('Average velocity')
    
    %plot average velocity as a function of the geometry
    figure
    plot(xi,Vmax,'o-')
    grid on
    xlabel('\xi')
    ylabel('V_{max}')
    title('Max velocity')
    
    %plot average velocity as a function of the geometry
    figure
    plot(xi,VmaxDim,'o-')
    grid on
    xlabel('\xi')
    ylabel('V_{max}[m/s]')
    title('Max velocity')

end

if CompareVisc==1

    %allocate cell size
    cellLegend = cell(size(visc));
    
    %loop in visc
    figure
    for k = 1:numel(visc)
        
        display([num2str(k) 'of ' num2str(numel(visc))])

        load(['ConicalMotor_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc(k)) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(10) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %store element in cell for legend
        cellLegend(k) = {['\lambda=' num2str(visc(k))]};
        
        %velocity data
        ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        for i = 1:ite
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end);

        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Micromotor velocity Ca=' num2str(Ca)])
        hold on

    end

    legend(cellLegend,'Location','Best')

end

if CompareHrep==1

    %allocate cell size
    cellLegend = cell(size(Hrep));
    
    %loop in visc
    figure
    for k = 1:numel(Hrep)
        
        display([num2str(k) 'of ' num2str(numel(Hrep))])

        load(['ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep(k)) '_IntenRep=' num2str(IntensRep) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(10) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %store element in cell for legend
        cellLegend(k) = {['H_{rep}=' num2str(Hrep(k))]};
        
        %velocity data
        ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        ForceDrop = zeros(1,ite);
        for i = 1:ite
            
            %useful variables
            m = find(risa(2:end,i)==risa(1,i));
            %indVel = find(risy(:,i)==0,2,'first');
            indXY = find(risb(:,i)==0,2,'first');
            aDrop = risa(m+2:indXY(2)-1,i);
            bDrop = risb(m+2:indXY(2)-1,i);
            
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
            
            %check force free condition
            ForceDrop(i) = forceOnDrop(aDrop,bDrop,KKK(1:nbr_el(i)+1,i),1);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end);
    
        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Motor vel Ca=' num2str(Ca) 'Resize dfX=' num2str(Resize)])
        hold on
        
        figure(2)
        plot(ttt,ForceDrop)
        grid on
        xlabel('t')
        ylabel('F')
        title(['Force on drop Ca=' num2str(Ca) 'Resize dfX=' num2str(Resize)])
        hold on

    end

    legend(cellLegend,'Location','Best')

end

if CompareResize==1

    %allocate cell size
    cellLegend = cell(size(Hrep));
    
    %loop in visc
    figure
    for k = 1:numel(Resize)
        
        display([num2str(k) 'of ' num2str(numel(Resize))])

        load(['ConicalMotor_Resize=' num2str(Resize(k)) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(10) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat']);

        %store element in cell for legend
        cellLegend(k) = {['Resize=' num2str(Resize(k))]};
        
        %velocity data
        ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        ForceDrop = zeros(1,ite);
        for i = 1:ite
            
            %useful variables
            m = find(risa(2:end,i)==risa(1,i));
            indVel = find(risy(:,i)==0,2,'first');
            indXY = find(risb(:,i)==0,2,'first');
            aDrop = risa(m+2:indXY(2)-1,i);
            bDrop = risb(m+2:indXY(2)-1,i);
            
            %indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
            
            %check force free condition
            ForceDrop(i) = forceOnDrop(aDrop,bDrop,KKK(1:nbr_el(i)+1,i),1);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end);
    
        figure(1)
        plot(ttt,v)
        grid on
        xlabel('t')
        ylabel('v')
        title(['Micromotor velocity Ca=' num2str(Ca)])
        hold on
        
        figure(2)
        plot(ttt,ForceDrop)
        grid on
        xlabel('t')
        ylabel('F')
        title(['Force on drop Ca=' num2str(Ca)])
        hold on

    end

    figure(1)
    legend(cellLegend,'Location','Best')
    
    figure(2)
    legend(cellLegend,'Location','Best')

end

if CompareIntensRep==1

    %allocate cell size
    cellLegend = cell(size(IntensRep));
    
    %loop in visc
    figure
    for k = 1:numel(IntensRep)
        
        display([num2str(k) 'of ' num2str(numel(IntensRep))])
        
        %clear previos data
        clear risa
        clear risb
        clear risy

        if isempty(deflationWall)
            if Hrep==5||Hrep==6
                filename = [source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep(k)) '_RepLenght=' num2str(repLenght(k)) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem(k)) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(10) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat'];
            else
                filename = [source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep(k)) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(10) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat'];
            end
        else
            filename = [source 'ConicalMotor_deflationWall=' num2str(deflationWall) '_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(Hrep) '_IntenRep=' num2str(IntensRep(k)) '_RepLenght=' num2str(repLenght(k)) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(elem(1)) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=' num2str(1) '_L=' num2str(Lenght) '_alpha=' num2str(alpha) '_InPos=' num2str(InPos) '_RK=' num2str(RK) '.mat'];
        end
            
        load(filename);
        
        %store element in cell for legend
        if Hrep==5||Hrep==6
            cellLegend(k) = {['A=' num2str(IntensRep(k)) ' \delta=' num2str(repLenght(k)) ' el=' num2str(elem(k))]};
        else
            cellLegend(k) = {['A=' num2str(IntensRep(k))]};
        end
        %velocity data
        ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1);
        v = zeros(1,ite);
        for i = 1:ite
            indVel = find(risy(:,i)==0,2,'first');
            v(i) = risy(indVel(2)-1,i);
        end
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end);

        plot(ttt,v)
        axis([0 time(end) 0 maxVel])
        grid on
        xlabel('t')
        ylabel('v')
        title(['Motor velocity Ca=' num2str(Ca) ' \lambda=' num2str(visc)])
        hold on

    end

    legend(cellLegend,'Location','Best')

end



