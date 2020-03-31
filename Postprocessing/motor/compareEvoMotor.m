%post processing micromotor

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%parameters
MotorFree = 0*ones(7,1);
Resize = 0*ones(8,1);
manyCa = 0.01*ones(8,1);
inflate = 0*ones(1,16);
dt = 0.01*ones(1,16);
theta = -0.05*ones(1,16);
alpha = 0.8*ones(1,16);
mass = 1*ones(1,16);
elem = 308*ones(1,16);
visc = 0*ones(1,16);
InPos = 1*ones(1,16);
Lenght = 10*ones(1,16);
Hrep = 6*ones(1,16);
IntensRep = -[100 10 1 0.1 100 10 1];
RepLenght = [0.1*ones(4,1); 0.01*ones(3,1)];
RK = 2*ones(1,16);
deflationWall = 1;

color = 'rgbkyrgbkyrgbkyrgbky';

%options
stop = 999;        %stop at this loop
step = 1;           %
plotRealTIME = 0;   %plot swimmer evolution in real time
expCompare = 2/3;   %compare with this scaling
stopCa = Inf;
res = 0;

%initialize
cellLegend = cell(numel(MotorFree));

for k = 1:numel(MotorFree)
    
%results path
here = pwd;
if res==1
    source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/';
elseif res==0
    source = '~/Documents/MATLAB/droplet_simulations/server/';
end
cd(source)

%load results
if isempty(deflationWall)
    filename = ['ConicalMotor_Resize=' num2str(Resize(k)) '_MotorFree=' num2str(MotorFree(k)) '_Hrep=' num2str(Hrep(k)) '_IntenRep=' num2str(IntensRep(k)) '_RepLenght=' num2str(RepLenght(k)) '_Inflate=' num2str(inflate(k)) '_theta=' num2str(theta(k)) '_el=' num2str(elem(k)) '_dt=' num2str(dt(k)) '_visc=' num2str(visc(k)) '_Ca=' num2str(manyCa(k)) '_R=' num2str(1) '_L=' num2str(Lenght(k)) '_alpha=' num2str(alpha(k)) '_InPos=' num2str(InPos(k)) '_RK=' num2str(RK(k)) '.mat'];
else
    filename = ['ConicalMotor_deflationWall=' num2str(deflationWall) '_Resize=' num2str(Resize(k)) '_MotorFree=' num2str(MotorFree(k)) '_Hrep=' num2str(Hrep(k)) '_IntenRep=' num2str(IntensRep(k)) '_RepLenght=' num2str(RepLenght(k)) '_Inflate=' num2str(inflate(k)) '_theta=' num2str(theta(k)) '_el=' num2str(elem(k)) '_dt=' num2str(dt(k)) '_visc=' num2str(visc(k)) '_Ca=' num2str(manyCa(k)) '_R=' num2str(1) '_L=' num2str(Lenght(k)) '_alpha=' num2str(alpha(k)) '_InPos=' num2str(InPos(k)) '_RK=' num2str(RK(k)) '.mat'];
end
breakNow = 0;

%clear data from previous plotting
clear risa
clear risb
clear risy

try
    load(filename)
catch
    warning('No data')
    breakNow = 1;
end
cd(here)

%store element in cell for legend
%cellLegend(k) = {['Compensation=' num2str(PARAM.resizeDFX) ' =' num2str(PARAM.q)]};

q = PARAM.q;
[~,loop] = size(risa);
vel = zeros(1,loop);
mesh = zeros(1,loop);
V = zeros(1,loop);
Xdrop = zeros(1,loop);
lubricationFilm = zeros(1,loop);
Ca = zeros(1,loop);
A = zeros(1,loop);
forceX = zeros(1,loop);
comp = zeros(1,loop);

%volume at first iteration
m = find(risa(2:end,1)==risa(1,1));
indXY = find(risb(:,1)==0,2,'first');
a_drop = risa(m+2:indXY(2)-1,1);
b_drop = risb(m+2:indXY(2)-1,1);
Vin = axis_int_gauss_vect(a_drop',b_drop');

broken = 0; %if the loop breaks or not
for i = 1:step:loop
    
    if breakNow==1
        break;
    end
    
    disp([num2str(i) ' of ' num2str(loop)])
    
    try
    m = find(risa(2:end,i)==risa(1,i));
    if PARAM.MotorFree==2
        indVel = find(risy(:,i)==0,3,'first');
    else
        indVel = find(risy(:,i)==0,2,'first');
    end
    indXY = find(risb(:,i)==0,2,'first');
    if PARAM.MotorFree==2
        vel(i) = risy(indVel(3)-1,i);
    else
        vel(i) = risy(indVel(2)-1,i);
    end
    mesh(i) = indXY(2)-1;
    a_wall = risa(1:m+1,i);
    b_wall = risb(1:m+1,i);
    
    if PARAM.MotorFree==2
        uDrop = risy(2*m+1:2:indVel(3)-3,i);
        vDrop = risy(2*m+2:2:indVel(3)-2,i);
    else
        uDrop = risy(2*m+1:2:indVel(2)-3,i);
        vDrop = risy(2*m+2:2:indVel(2)-2,i);
    end
    
    a_drop = risa(m+2:indXY(2)-1,i);
    b_drop = risb(m+2:indXY(2)-1,i);
    
    if i==1
        aWallFirst = a_wall;
        bWallFirst = b_wall;
        aDropFirst = a_drop;
        bDropFirst = b_drop;
    elseif i==loop||i==stop
        aWallLast = a_wall;
        bWallLast = b_wall;
        aDropLast = a_drop;
        bDropLast = b_drop;
    end
    
    qNow = numel(a_drop)-1;
    
    phi = FigInOut(a_drop',b_drop',a_wall,b_wall);
    comp(i) = max(phi>pi);
    try
        V(i) = axis_int_gauss_vect(a_drop',b_drop');
    catch
        warning('No data')
        broken=1;
        breakLoop = i-step;
        break;
    end
    
    A(i) = surf_gauss_vect(a_drop',b_drop');
    Xdrop(i) = center_mass(a_drop',b_drop');
    Un = DropNormalVelocity(a_drop',b_drop',uDrop,vDrop);
    Ca(i) = abs(DropVelocityAxis(a_drop',b_drop',Un));
    lubricationFilm(i) = computeDistance(a_drop(round(numel(a_drop)/2)),b_drop(round(numel(a_drop)/2)),a_wall,b_wall);
    
    %consider if the motor is moving
    if PARAM.MotorFree~=2
        Xdrop(i) = Xdrop(i)-max(a_wall);
        Ca(i) = Ca(i)-vel(i);
    end
    
    if plotRealTIME==1
    
    figure(1)
    subplot(2,2,1)
    plot(a_wall,b_wall,'kx-',a_wall,-b_wall,'k-')
    hold on
    plot(a_drop,b_drop,'rx-',a_drop,-b_drop,'r-')
    hold off
    axis equal
    grid on
    xlabel('x')
    ylabel('r')
    axis([-5 11 -2.5 2.5])
    %axis([-0.2 0.4 0.95 1.15])
    title(['Ca=' num2str(PARAM.Ca) ' \alpha=' num2str(PARAM.alpha) ' \lambda=' num2str(PARAM.visc)])
    
    %plot motor velocity
    %figure(2)
    subplot(2,2,2)
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,vel(i),'ob')
    drawnow
    grid on
    hold off
    xlabel('time')
    ylabel('velocity')
    
    %capillary numner
    subplot(2,2,3)
    %hold on
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,abs(Xdrop(i)-Xdrop(1)),'ob')
    grid on
    hold off
    %hold off
    xlabel('t')
    ylabel('x_{pos}')
    
    %drop position
    subplot(2,2,4)
    %hold on
    hold on
    plot(i*PARAM.checkpoint*PARAM.deltaT,Ca(i),'ob')
    grid on
    hold off
    %hold off
    xlabel('t')
    ylabel('Ca')
    
%     %drop volume
%     subplot(2,2,4)
%     %hold on
%     hold on
%     plot(i*PARAM.checkpoint*PARAM.deltaT,V(i),'ob')
%     grid on
%     hold off
%     %hold off
%     xlabel('x')
%     ylabel('A')
    
    end
    
    catch
        warning('No data')
        vel = vel(1:i);
        Ca = Ca(1:i);
        lubricationFilm = lubricationFilm(1:i);
        Xdrop = Xdrop(1:i);
        aWallLast = a_wall;
        bWallLast = b_wall;
        aDropLast = a_drop;
        bDropLast = b_drop;
        break;
    end
    
%     if i>stop
%         disp('STOP ARTIFICIALLY!')
%         vel = vel(1:i);
%         Ca = Ca(1:i);
%         lubricationFilm = lubricationFilm(1:i);
%         Xdrop = Xdrop(1:i);
%         break;
%     end

    
end

%plot motor velocity
if PARAM.MotorFree~=2
    figure(2)
    if k~=1
        hold on
    end
    plot((0:numel(Xdrop)-1)*PARAM.deltaT*PARAM.checkpoint,vel,'-')
    drawnow
    grid on
    hold off
    xlabel('time')
    ylabel('velocity')
end

%plot first and last motor position
figure(3)
subplot(2,1,1)
if k~=1
        hold on
end
plot(aWallFirst,bWallFirst,'k')
hold on
plot(aWallFirst,-bWallFirst,'k')
plot(aDropFirst,bDropFirst,color(k))
plot(aDropFirst,-bDropFirst,color(k))
axis equal
axis([-5 PARAM.L+1 -2.5 2.5])
xlabel('x')
ylabel('r')
grid on
title('Initial shape')
subplot(2,1,2)
plot(aWallLast,bWallLast,'k')
hold on
plot(aWallLast,-bWallLast,'k')
plot(aDropLast,bDropLast,color(k))
plot(aDropLast,-bDropLast,color(k))
xlabel('x')
ylabel('r')
grid on
axis equal
axis([-5 PARAM.L+1 -2.5 2.5])
title('Final shape')

%drop position
figure(4)
if k~=1
        hold on
end
%subplot(1,3,2)
loglog((0:numel(Xdrop)-1)*PARAM.deltaT*PARAM.checkpoint,abs(Xdrop-Xdrop(1)),'x')
grid on
hold off
%hold off
xlabel('t')
ylabel('x_{pos}')
title('Drop position')
    
%capillary number
figure(5)
if k~=1
        hold on
end
loglog((0:numel(Ca)-1)*PARAM.deltaT*PARAM.checkpoint,Ca,'x')
grid on
hold off
%hold off
xlabel('t')
ylabel('Ca')
title('Capillary number')

figure(6)
%subplot(1,3,3)
if k~=1
        hold on
end
if stopCa~=Inf
    ind = find(Ca>6e-3,1,'first');
else
    ind = numel(Ca);
end
loglog(Ca(1:ind),lubricationFilm(1:ind),'x')
hold on
if k==numel(MotorFree)
loglog(Ca(1:ind),Ca(1:ind).^expCompare,'k')
end
grid on
hold off
xlabel('Ca')
%leg = ['Ca^{' num2str(expCompare) '}'];
%legend('simulation',leg)
ylabel('h')
title(['Lubrication film compared to Ca^{' num2str(expCompare) '}'])

end

