%compare quantities from the simulation and from the model

clear variables
close all

%problem parameters
%source = '~/Documents/MATLAB/droplet_simulations/results/micromotor/';
source = '~/Documents/MATLAB/droplet_simulations/server/';
dest = '~/Documents/research_notes/EuroMecSevilla/presentation/movies/frames/';
theta = -0.1;
Ca = 0.0;
dt = 0.01;
visc = 0.1;
inflate = 0;
element = 308;
MotorFree = 0;
TypeRep = 6;
IntensRep = -1e-1;
RepLenght = 1e-1;
Resize = 0;
alpha = 1.2;
L = 10;
Inpos = 5;

%upload data
filename = [source 'ConicalMotor_Resize=' num2str(Resize) '_MotorFree=' num2str(MotorFree) '_Hrep=' num2str(TypeRep) '_IntenRep=' num2str(IntensRep) '_RepLenght=' num2str(RepLenght) '_Inflate=' num2str(inflate) '_theta=' num2str(theta) '_el=' num2str(element) '_dt=' num2str(dt) '_visc=' num2str(visc) '_Ca=' num2str(Ca) '_R=1_L=' num2str(L) '_alpha=' num2str(alpha) '_InPos=' num2str(Inpos) '_RK=2.mat'];
load(filename)

%plotting options
step = 100;                                                      %plot every few iteration
SubPlot = 0;                                                    % plot with subplot
ComputeDropVel = 1;                                             % compute velocity of the droplet from the solution
ComputeDropPos = 1;                                             % compute droplet position
ComputeDropSurface = 1;                                         % compute surface area of the bubble
ComputeDropVolume = 1;                                         % compute volume of the bubble
PlotDropVel = 0;                                                % plot velocity of the droplet from the solution
PlotVelNow = 0;                                                 % plot motor velocity at this instant
computeFilmThickness = 1;

%till this time
time = 5800;

%first coordinate of motor
m = find(risa(2:end,1)==risa(1,1));
indNode = find(risb(:,1)==0,2,'first');
aMotorFirst = risa(1:m+1,1);
bMotorFirst = risb(1:m+1,1);

%initialize
ite = round(time(end)/PARAM.checkpoint/PARAM.deltaT+1)-1;
v = zeros(1,ite/step);
displace = zeros(1,ite/step);
vDrop = zeros(1,ite/step);
dropPos = zeros(1,ite/step);
AreaDrop = zeros(1,ite/step);
xcmDrop = zeros(1,ite/step);
volumeDrop = zeros(1,ite/step);
AreaModel = zeros(1,ite/step);
filmThick = zeros(1,ite/step);

for k = 1:ite/step
    %display(i)
    i = k*step+1;
    indVel = find(risy(:,i)==0,2,'first');
    v(k) = risy(indVel(2)-1,i);
    displace(k) = risa(1,i)-risa(1,1);
    
    if ComputeDropVel==1||ComputeDropPos==1
        
        display([num2str(i) ' of ' num2str(ite)])
        
        %often used
        m = find(risa(2:end,i)==risa(1,i));
        indNode = find(risb(:,i)==0,2,'first');
        q = indNode(2)-m-3;
        aDrop = risa(m+2:m+q+2,i);    bDrop = risb(m+2:m+q+2,i);
        xWall = risa(1:m+1,i);    yWall = risb(1:m+1,i);
        
        if ComputeDropVel==1
            VelInterfaceX = risy(2*m+1:2:2*(m+q)+1,i);
            VelInterfaceY = risy(2*m+2:2:2*(m+q)+2,i);

            %compute vector normal to the interface
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(aDrop',bDrop');
            N = [by./sqrt(bx.*bx+by.*by) (by(end)+2*cy(end)+3*dy(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)));...
                -bx./sqrt(bx.*bx+by.*by) (-bx(end)-2*cx(end)-3*dx(end))/sqrt((bx(end)+2*cx(end)+3*dx(end))*(bx(end)+2*cx(end)+3*dx(end))+(by(end)+2*cy(end)+3*dy(end))*(by(end)+2*cy(end)+3*dy(end)))];

            %velocity normal to interface
            Unormal = VelInterfaceX'.*N(1,:) + VelInterfaceY'.*N(2,:);
            
            %compute velocity
            vDrop(k) = DropVelocityAxis(aDrop',bDrop',Unormal');
        
        end
        
        if ComputeDropPos==1
            
            dropPos(k) = center_mass(aDrop,bDrop);
            
        end
        
        if ComputeDropSurface==1
            
            %compute surface area
            AreaDrop(k) = surf_gauss_vect(aDrop',bDrop');
            %xcmDrop(k) = center_mass_gauss(aDrop',bDrop');
        
        end
        
        if ComputeDropVolume == 1
            
            %compute volume
            volumeDrop(k) = axis_int_gauss_vect(aDrop',bDrop');
            
        end
        
        if computeFilmThickness==1
            
            x0 = aDrop(round(numel(aDrop)/2));
            y0 = aDrop(round(numel(bDrop)/2));
            filmThick(k) = min(sqrt((xWall-x0).^2+(yWall-y0).^2));
            
        end
        
    end
    
    %area from the model with better geometry
    PARAM.h = 0.15;
    [AreaModel(k),R1,R2,H] = findR1R2HforModel(PARAM.h,PARAM.L,PARAM.R,volumeDrop(k),dropPos(k),PARAM.theta);
    
    %plot current model shape
    [xDrop,yDrop] = draw_test3(dropPos(k),0,H,R1,R2,100,PARAM.theta);
    
    %checkArea(k) = surf_gauss_vect(xDrop,yDrop);
    
    %plot drop
    figure(1)
    plot([xDrop flip(xDrop)],[yDrop -flip(yDrop)])
    hold on
    plot(xWall,yWall,'k')
    plot(xWall,-yWall,'k')
    xlabel('x')
    ylabel('r')
    grid on
    axis equal
    axis([-0.5 10.5 -4 4])
    hold off
    drawnow
    
end
tttHere = 0:PARAM.checkpoint*PARAM.deltaT*step:time(end)-PARAM.checkpoint*PARAM.deltaT;

%compute area and area variation of the droplet as in the model
VolumeModel = 4/3*pi*PARAM.alpha^3;
Li = -(1-sin(PARAM.theta))/tan(PARAM.theta);
%AreaModel = -2*VolumeModel/PARAM.theta./(10-dropPos+Li);
%dAreaDxModel = -2*VolumeModel/PARAM.theta./(10-dropPos+Li).^2;
dAreaDxModel = -diff(AreaModel)./diff(10-dropPos);

%motor friction
Fdrag = (2*pi*L*v)/(log(L/2/radius)+log(2)-0.5);
figure
plot(tttHere,Fdrag,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('Drag on motor')

figure
plot(tttHere,v,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('v_{motor}')

figure
plot(tttHere,dropPos,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('pos_{bubble}')
    
figure
plot(10-dropPos(1:end-1),diff(AreaDrop)./diff(dropPos),'LineWidth',2)
hold on
plot(10-dropPos(1:end-1),dAreaDxModel,'LineWidth',2)
legend('Simulations','Model','Location','Best')
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('space')
ylabel('dA_{bubble}/dx')

figure
plot(tttHere,vDrop,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('V_{bubble}')

figure
errV = (volumeDrop-volumeDrop(1))/volumeDrop(1);
plot(tttHere,errV,'LineWidth',2)
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('errV')

figure
plot(tttHere,AreaDrop,'LineWidth',2)
hold on
plot(tttHere,AreaModel,'LineWidth',2)
legend('Simulations','Model','Location','Best')
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('A_{bubble}')
    
figure
plot(tttHere(1:end-1),diff(AreaDrop)./diff(tttHere),'LineWidth',2)
hold on
plot(tttHere(1:end-1),diff(AreaModel)./diff(tttHere),'LineWidth',2)
legend('Simulations','Model','Location','Best')
title(['Ca = ' num2str(PARAM.Ca) ' \lambda = ' num2str(PARAM.visc) ' \theta=' num2str(PARAM.theta)])
hold on
grid on
xlabel('time')
ylabel('dA_{bubble}/dt')

figure
loglog(abs(vDrop),filmThick)
hold on
loglog(abs(vDrop),abs(vDrop).^(2/3)+filmThick(end),'k')
xlabel('Ca')
ylabel('h')
grid on
title('Film thickness')











