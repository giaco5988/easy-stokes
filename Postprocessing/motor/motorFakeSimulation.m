%plot DNS at some desired times

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

close all
clear variables

%problem parameters
source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/';
here = pwd;
BEM = '~/Documents/MATLAB/droplet_simulations/conicalMotor';
%source = '~/Documents/MATLAB/droplet_simulations/server/';
dest = '~/Documents/research_notes/EuroMecSevilla/presentation/movies/frames/';
theta = -0.02;
Ca = 0.01;
dt = 0.01;
visc = 0.1;
inflate = 0;
element = 308;
MotorFree = 0;
TypeRep = 6;
IntensRep = -1e-1;
RepLenght = 1e-1;
Resize = 0;
alpha = 0.8;
L = 10;
Inpos = L/2;

%upload data
filename = [source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(TypeRep) '_IntenRep=' num2str(IntensRep) '_RepLenght=' num2str(RepLenght) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(element) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=' num2str(L) '_alpha=' num2str(alpha) '_InPos=' num2str(Inpos) '_RK=2.mat'];
load(filename)

%define new parameters
PARAMhere.Ca = Ca;
PARAMhere.theta = theta;
PARAMhere.visc = visc;
PARAMhere.inflate = inflate;
PARAMhere.element = element;
PARAMhere.MotorFree = MotorFree;
PARAMhere.repulsive = 0;
PARAMhere.coeffRepulsive = IntensRep;
PARAMhere.repulsiveOn = RepLenght;
PARAMhere.resizeDFX = Resize;
PARAMhere.visc2 = 1;
PARAMhere.cfunction = PARAM.cfunction;
PARAMhere.spline = PARAM.spline;
PARAMhere.doubleLayerSingTreat = PARAM.doubleLayerSingTreat;
PARAMhere.NS_wall_drop = PARAM.NS_wall_drop;
PARAMhere.NS_drop_wall = PARAM.NS_drop_wall;

%plotting options
ComputeDropVel = 1;                                             % compute velocity of the droplet from the solution
semilog = 0;

%which time to plot
time = linspace(10,5000,50);

%initialization
vDrop = zeros(numel(time),1);
vDropHere = zeros(numel(time),1);
vMotor = zeros(numel(time),1);
vMotorHere = zeros(numel(time),1);

for i = 1:numel(time)

        display([num2str(i) ' of ' num2str(numel(time))])

        if i==1
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1);
        else
            ite = round(time(i)/PARAM.checkpoint/PARAM.deltaT+1)-1;
        end

        %often used
        m = find(risa(2:end,ite)==risa(1,ite));
        indNode = find(risb(:,ite)==0,2,'first');
        q = indNode(2)-m-3;
        
        %velocity from original simulation
        indVel = find(risy(:,ite)==0,2,'first');
        vMotor(i) = risy(indVel(2)-1,ite);

        %coordinates
        aMotor = risa(1:m+1,ite);       bMotor = risb(1:m+1,ite);
        aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);
        a = [aMotor; aDrop]';     b = [bMotor; bDrop]';
        
        %compute original velocity
        if ComputeDropVel==1
    
            % velocities on the interface
            VelInterfaceX = risy(2*m+1:2:2*(m+q)+1,ite);
            VelInterfaceY = risy(2*m+2:2:2*(m+q)+2,ite);

            %compute vector normal to the interface
            aDrop = risa(m+2:m+q+2,ite);    bDrop = risb(m+2:m+q+2,ite);
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(aDrop',bDrop');
            N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
                -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

            %velocity normal to interface
            Unormal = VelInterfaceX'.*N(1,:) + VelInterfaceY'.*N(2,:);

            %compute velocity
            vDrop(i) = DropVelocityAxis(aDrop',bDrop',Unormal');
        
        end
        
        %run simaultion with modified parameters
        PARAMhere.Xinj = center_mass(aDrop,aDrop);
        PARAMhere.m = m;    PARAMhere.q = q;
        PARAMhere.Qsource = PARAMhere.Ca;
        cd(BEM)
        [yhere,K,K1,K2,N] = BEMmotor(a,b,PARAMhere);
        cd(here)
        
        vMotorHere(i) = yhere(end);
        
        if ComputeDropVel==1
        
            %often used
            VelInterfaceX = yhere(2*m+1:2:2*(m+q)+1);
            VelInterfaceY = yhere(2*m+2:2:2*(m+q)+2);

            %compute vector normal to the interface
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(aDrop',bDrop');
            N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
                -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

            %velocity normal to interface
            Unormal = VelInterfaceX'.*N(1,:) + VelInterfaceY'.*N(2,:);

            %compute velocity
            vDropHere(i) = DropVelocityAxis(aDrop',bDrop',Unormal');
        
        end
            
end

figure
if semilog==0
    plot(time,vMotor)
elseif semilog==1
    semilogy(time,vMotor)
elseif semilog==2
    loglog(time,vMotor)
end
hold on
if semilog==0
    plot(time,vMotorHere,'k')
elseif semilog==1
    semilogy(time,vMotorHere,'k')
elseif semilog==2
    loglog(time,vMotorHere,'k')
end
grid on
xlabel('t')
ylabel('V_{motor}')
title(['Motor velocity \theta=' num2str(theta) ' Ca=' num2str(Ca)])
legend('Original simulation','Re-run simulation','Location','Best')

figure
if semilog==0
    plot(time,vDrop)
elseif semilog==1
    semilogy(time,vDrop)
elseif semilog==2
    loglog(time,vDrop)
end
hold on
if semilog==0
    plot(time,vDropHere,'k')
elseif semilog==1
    semilogy(time,vDropHere,'k')
elseif semilog==2
    loglog(time,vDropHere,'k')
end
grid on
xlabel('t')
ylabel('V_{drop}')
title(['Drop velocity \theta=' num2str(theta) ' Ca=' num2str(Ca)])
legend('Original simulation','Re-run simulation','Location','Best')


















