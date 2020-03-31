%plot DNS at some desired times

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

close all
clear variables

%problem parameters
source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/parametric_study/theta005/';
dest = '~/Documents/research_notes/someMovies/motor/frames/';
theta = -0.05;
Ca = 0.001;

%upload data
filename = [source 'LabFrame_mass=1_theta=' num2str(theta) '_el=320_dt=0.0001_visc=1_Ca=' num2str(Ca) '_R=1_L=10_alpha=0.8_InPos=4_RK2.mat'];
load(filename)

%plotting options
PlotSnapshot = 1;       % plot snapshots
SubPlot = 1;            % plot with subplot
PlotVelField = 0;       % plot velocity field
SaveVelField = 0;       % save vel field
FrameLab = 1;           % 0 plot in motor frame, 1 lab frame or 2 bubble frame

%geometry option
cutX = 8;   cutY = 3;   shift = 4;
xRange1 = -cutX + shift;   xRange2 = cutX + shift;
yRange1 = 0.0;   yRange2 = cutY;
MeshFine = 4;

%which time to plot
time = linspace(0,20,12);
%time = [0 7.5 16 18.5 20];
%time = [5 11 13.7 14.5];
%time = [0 11 13.7 14.5];
%time = [11 12];
time2 = [0 3 6 20]*1;

%when I have run two separate simulations
PASTE = 0;
startTime = 0;

if PASTE==0
    close all
end

if PASTE==1
    time=time2;
end

%figure in pixel
width = 1400;
height = 900;

%velocity data
ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1)-1;
v = zeros(1,ite);
for i = 1:ite
    %display(i)
    indVel = find(risy(:,i)==0,2,'first');
    v(i) = risy(indVel(2)-1,i);
end
ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(end)-PARAM.checkpoint*PARAM.deltaT;

if PASTE==0
    figure(2)
    plot(ttt,v,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
elseif PASTE==1
    figure(2)
    hold on
    plot(ttt+startTime,v,'LineWidth',2)
    title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
end
hold on
grid on
xlabel('time')
%axis([0 time(end) min(risy(m+q+2,:)) max(risy(m+q+2,:))])
ylabel('v')

if PlotSnapshot==1
    for i = 1:numel(time)

        display(i)

        if i==1
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1);
        else
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1)-1;
        end

        %often used
        m = PARAM.m;    %q = nbr_el(ite);
        indNode = find(risb(:,ite)==0,2,'first');
        q = indNode(2)-m-3;

        aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
        aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);

        %velocity data
        vnow = v(1:ite);
        ttt = 0:PARAM.checkpoint*PARAM.deltaT:time(i);

        if PASTE==0
            figure(1)
            fig = figure(1);
        elseif PASTE==1
            figure(5)
            fig = figure(5);
        end
        
        if SubPlot==1
            fig.Position = [400 200 width height];
            subplot(numel(time)/4,4,i)
        else
            figure(1)
            %hold off
        end
        
        if PlotVelField==1
            
            x = linspace(xRange1,xRange2,MeshFine*(xRange2-xRange1));
            y = linspace(yRange1,yRange2,MeshFine*yRange2);
            [X,Y] = meshgrid(x,y);
            
            [U,V,uFrame] = velocity_field_motor_tracers([aMotor; aDrop]',[bMotor; bDrop]',risy,X,Y,ite,PARAM,FrameLab,0,1);
            
            %U = U+uFrame;
%             subplot(2,1,1)
%             plot(time(1:i),v(1:i))
%             grid on
%             xlabel('t')
%             ylabel('V')
%             title('Motor Velocity')
            
            %figure(1)
            %subplot(2,1,2)
            %normalizeVel = sqrt((U-uFrame).^2+V.^2);
            %quiver(X,-Y,(U-uFrame)./normalizeVel,-V./normalizeVel,'b')
            quiver(X,-Y,(U-uFrame),-V,'b')
            hold on
            streamslice(X,Y,U-uFrame,V,'b')
            
        end
        plot(aDrop,bDrop,'r',aDrop,-bDrop,'r')
        hold on
        plot(aMotor,bMotor,'k',aMotor,-bMotor,'k')
        grid on
        axis equal
        axis([xRange1 xRange2 -yRange2 yRange2])
        xlabel('x')
        ylabel('r')
        drawnow
        hold off
        
        
        if PASTE==0
            
            if SubPlot==0
                if FrameLab==0
                    title(['Flow field motor frame, t=' num2str(time(i))])
                elseif FrameLab==1
                    title(['Flow field lab frame, t=' num2str(time(i))])
                elseif FrameLab==2
                    title(['Flow field droplet frame, t=' num2str(time(i))])
                end
            end
            
        elseif PASTE==1
            title(['t=' num2str(time(i)+startTime)])
        end
        
        if SaveVelField==1

            
            fig = gcf;
            fig.Position = [400 200 width height];
            name = 'MotorFrameVelFieldMask';
            print('-dpng','-loose','-r100',[dest name sprintf('%04d',i) '.png'])
            
            
        end

        
        if PASTE==0
            figure(2)
            hold on
            plot(ttt(end),vnow(end),'.','MarkerSize',40)
        elseif PASTE==1
            figure(4)
            hold on
            plot(ttt(end)+startTime,vnow(end),'.','MarkerSize',40)
        end

    end

    if PASTE==0
        
        if numel(time)==4
            legend('v',['t=' num2str(time(1))],['t=' num2str(time(2))],['t=' num2str(time(3))],['t=' num2str(time(4))],'Location','Best')
        elseif numel(time)==3
            legend('v',['t=' num2str(time(1))],['t=' num2str(time(2))],['t=' num2str(time(3))],'Location','Best')
        end
        
        title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta) ' x_0=' num2str(PARAM.start)])
    elseif PASTE==1
        legend('v',['t=' num2str(time(1)+startTime)],['t=' num2str(time(2)+startTime)],['t=' num2str(time(3)+startTime)],['t=' num2str(time(4)+startTime)],'Location','Best')
        title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta) ' x_0=' num2str(PARAM.start)])
    end

end