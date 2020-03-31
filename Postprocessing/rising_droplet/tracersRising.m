%compute moving tracers for velocity field

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%load data
load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/prolate/q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.07_Ca=6_RK2.mat')
dest = '~/Documents/research_notes/someMovies/relaxation/';

%range to plot
xRange = 2.5;
yRange = 2.5;
MeshFine = 35;
%how many loop
loop = 170;
step = 160;
DistReplace = Inf;
dropFrame = 1;

%options
plotOnlyShape = 0;
plotVelField = 1;
saveRadialVelField = 0;
saveAxialVelField = 0;
saveQuiverVelField = 0;
saveIntensitylVelField = 0;
plotTracers = 0;

%only inside the droplet cartesian
resX = 10;
resY = 10;
a = risa(:,1)';
b = risb(:,1)';
x = linspace(min(a)-0.1*min(a),max(a)-0.1*max(a),2*resX);
y = linspace(0,max(b),resY);
[X,Y] = meshgrid(x,y);
InOut = FigInOut(a,b,X(:),Y(:));
InOut = reshape(InOut,numel(y),numel(x));
X = (InOut>pi).*X + (InOut<pi)*0;
Y = (InOut>pi).*Y + (InOut<pi)*0;

%count effectives loops
iteSave = 0;

for i = 1:step:loop
    
    iteSave = iteSave+1;
    
    disp([num2str(i) ' of ' num2str(loop)])
    
    ite = i;
    
    a = risa(:,ite)';
    b = risb(:,ite)';
    
    if dropFrame==1
        xcm = center_mass(a,b);
        a = a-xcm;
    end
    
    if plotOnlyShape==1
        
        %plot shape
        figure(1)
        plot(b,-a,'k',-b,-a,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        grid on
        xlabel('x')
        ylabel('r')
        drawnow
        
    end
    
    if plotTracers==1
    
        %plot tracers
        figure(1)
        plot(Y,-X,'.b','MarkerSize',10)
        hold on
        plot(-Y,-X,'.b','MarkerSize',10)
        plot(b,-a,'k',-b,-a,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        grid on
        xlabel('x')
        ylabel('r')
        drawnow
        hold off

        %compute velovity
        [U,V] = velocity_field_buoyancy_tracers(risa(:,ite)',risb(:,ite)',risy(:,ite),visc,Ca,X,Y,DistReplace,dropFrame);

        %update tracers
        X = X+deltaT*U*step;
        Y = Y+deltaT*V*step;
    
    end 
        
    if plotVelField==1
        
        finess = MeshFine;
        x = linspace(-xRange,xRange,finess);
        y = linspace(0,yRange,finess);
        [X,Y] = meshgrid(x,y);
        
        %compute velocity
        [U,V] = velocity_field_buoyancy_tracers(risa(:,ite)',risb(:,ite)',risy(:,ite),visc,Ca,X,Y,DistReplace,dropFrame);
        
        %get colorbar limits
        if i==1
           
            Umax = max(max(U));
            Vmax = max(max(V));
            IntensMax = max(max(sqrt(U.^2+V.^2)));
            
            Umin = min(min(U));
            Vmin = min(min(V));
            IntensMin = min(min(sqrt(U.^2+V.^2)));
            
        end
        
        %plot vel field
        figure(2)
        %contourf(X,Y,V,10,'LineStyle','none')
        contourf(Y,-X,V,10)
        hold on
        %contourf(X,-Y,V,10,'LineStyle','none')
        contourf(-Y,-X,V,10)
        plot(b,-a,'k',-b,-a,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        %grid on
        xlabel('x')
        ylabel('r')
        hold off
        title('Radial velocity field')
        colorbar
        caxis([Vmin Vmax])
        drawnow
        
        
        if saveRadialVelField==1

            name = 'relaxationRadialVelField';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
            
            
        end
        
        finess = MeshFine;
        x = linspace(-xRange,xRange,finess);
        y = linspace(0,yRange,finess);
        [X,Y] = meshgrid(x,y);
        
        %compute velocity
        %[U,V] = velocity_field_relaxation_tracers(risa(:,ite)',risb(:,ite)',risy(:,ite),visc,Ca,X,Y);
        
        %plot vel field
        figure(3)
        contourf(Y,-X,-U,10,'LineStyle','none')
        hold on
        contourf(-Y,-X,-U,10,'LineStyle','none')
        plot(a,b,'k',a,-b,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        %axis off
        %grid on
        xlabel('x')
        ylabel('r')
        hold off
        title('Axial velocity field')
        colorbar
        caxis([Umin Umax])
        drawnow
        
        
        if saveAxialVelField==1

            name = 'relaxationAxialVelField';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
            
            
        end
        
        %plot vel field
        figure(5)
        contourf(Y,-X,sqrt(U.^2+V.^2),10,'LineStyle','none')
        %contourf(X,Y,sqrt(U.^2+V.^2),10)
        hold on
        contourf(-Y,-X,sqrt(U.^2+V.^2),10,'LineStyle','none')
        %streamslice(X,-Y,U,-V,'b')
    
        plot(a,b,'k',a,-b,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        %grid on
        xlabel('x')
        ylabel('r')
        hold off
        title('Velocity field intensity')
        colorbar
        caxis([IntensMin IntensMax])
        drawnow
        
        %plot vel field
        figure(6)
        plot(b,-a,'k',-b,-a,'k')
        streamslice(Y',-X',V',-U')
        %hold on
        streamslice(-Y',-X',-V',-U')
        axis equal
        axis([-xRange xRange -yRange yRange])
        xlabel('x')
        ylabel('r')
        title('Velocity field intensity')
        drawnow
        hold off
        
        
        if saveIntensitylVelField==1

            name = 'relaxationIntensityVelField';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
            
            
        end
        
    end
    
end