%compute moving tracers for velocity field

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%load data
load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/relaxation/Relaxation_q=200_visc=1_Dt=0.01_loop=3000_DELTA=1.5_Ca=1_RK2.mat')
dest = '~/Documents/research_notes/someMovies/relaxation/';

%range to plot
xRange = 3.5;
yRange = 3.5;
MeshFine = 50;
%how many loop
loop = 10;
step = 10;

%options
plotOnlyShape = 1;
plotVelField = 0;
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
    
    if plotOnlyShape==1
        
        %plot shape
        figure(1)
        plot(a,b,'k',a,-b,'k')
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
        plot(X,-Y,'.b','MarkerSize',10)
        hold on
        plot(X,Y,'.b','MarkerSize',10)
        plot(a,b,'k',a,-b,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        grid on
        xlabel('x')
        ylabel('r')
        drawnow
        hold off

        %compute velovity
        [U,V] = velocity_field_relaxation_tracers(risa(:,ite)',risb(:,ite)',risy(:,ite),visc,Ca,X,Y);

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
        [U,V] = velocity_field_relaxation_tracers(risa(:,ite)',risb(:,ite)',risy(:,ite),visc,Ca,X,Y);
        
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
        contourf(X,Y,V,10)
        hold on
        %contourf(X,-Y,V,10,'LineStyle','none')
        contourf(X,-Y,V,10)
        plot(a,b,'k',a,-b,'k')
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
        contourf(X,Y,U,10,'LineStyle','none')
        hold on
        contourf(X,-Y,U,10,'LineStyle','none')
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
%         figure(4)
%         quiver(X,Y,U,V)
%         hold on
%         quiver(X,-Y,U,-V)
%         plot(a,b,'k',a,-b,'k')
%         axis equal
%         axis([-xRange xRange -yRange yRange])
%         grid on
%         xlabel('x')
%         ylabel('r')
%         hold off
%         title('Velocity field')
%         drawnow
%         
%         
%         if saveQuiverVelField==1
% 
%             name = 'relaxationQuiverVelField';
%             print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
%             
%             
%         end
        
        %plot vel field
        figure(5)
        contourf(X,Y,sqrt(U.^2+V.^2),10,'LineStyle','none')
        %contourf(X,Y,sqrt(U.^2+V.^2),10)
        hold on
        contourf(X,-Y,sqrt(U.^2+V.^2),10,'LineStyle','none')
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
        
        
        if saveIntensitylVelField==1

            name = 'relaxationIntensityVelField';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
            
            
        end
        
    end
    
end