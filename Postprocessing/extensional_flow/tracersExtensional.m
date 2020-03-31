%compute moving tracers for velocity field

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

%load data
load('/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/extensional/ForMovies/Extensional_q=100_visc=1_Dt=0.05_loop=10000_Ca=0.1_VR=1_RK2.mat')
dest = '~/Documents/research_notes/someMovies/extensional/';

%range to plot
xRange = 3.5;
yRange = 3.5;
MeshFine = 50;
%how many loop
LoopStart = 1;
LoopFinish = 2000;
step = 1;

%options
plotInterface = 0;
plotVelField = 1;
saveRadialVelField = 0;
saveAxialVelField = 0;
saveQuiverVelField = 0;
saveIntensitylVelField = 0;
plotTracers = 0;
plotVelFieldTrace = 0;
saveTracers = 0;
PlotSnapshots = 1;  Tsnap = 10;

%only inside the droplet cartesian
resX = 5;
resY = 40;
a = risa(:,1)';
b = risb(:,1)';
x = linspace(min(a)-0.1*min(a),max(a)-0.1*max(a),2*resX);
y = linspace(0.05,max(b),resY);
[X,Y] = meshgrid(x,y);
InOut = FigInOut(a,b,X(:),Y(:));
InOut = reshape(InOut,numel(y),numel(x));
X = (InOut>pi).*X + (InOut<pi)*0;
Y = (InOut>pi).*Y + (InOut<pi)*0;
% finess = MeshFine;
% x = linspace(-xRange,xRange,finess);
% y = linspace(0,yRange,finess);
% [X,Y] = meshgrid(x,y);

%count effectives loops
iteSave = 0;

%counters
CountSnap = 1;

for i = LoopStart:step:LoopFinish
    
    iteSave = iteSave+1;
    
    disp([num2str(i) ' of ' num2str(loop)])
    
    ite = i;
    
    a = risa(:,ite)';
    b = risb(:,ite)';
    
    %plot interface
    if plotInterface==1
        
        %compute the spline coeff
        [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(a,b);

        %compute splines coordinates
        t = 0:0.1:0.9;

        %INTERFACE
        ttt = repmat(t,1,numel(ax));
        axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
        bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
        cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
        dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
        ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
        byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
        cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
        dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));

        %splines coordinates
        xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a(end)];
        yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b(end)];
        
        figure(6)
        plot(xxx,yyy,'k')
        hold on
        plot(xxx,-yyy,'k')
        hold off
        xlabel('x')
        ylabel('r')
        axis equal
        grid on
        axis([-xRange xRange -yRange yRange])
        drawnow
        
    end
    
    if PlotSnapshots==1;
        
        if CountSnap>numel(Tsnap)
            break;
        end
            
        if (i-1)*deltaT==Tsnap(CountSnap)
            
            %compute the spline coeff
            [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric(a,b);

            %compute splines coordinates
            t = 0:0.1:0.9;

            %INTERFACE
            ttt = repmat(t,1,numel(ax));
            axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
            bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
            cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
            dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
            ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
            byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
            cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
            dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));

            %splines coordinates
            xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 a(end)];
            yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 b(end)];
            
            figure(7)
            subplot(numel(Tsnap),1,CountSnap)
            plot(xxx,yyy,'k')
            hold on
            plot(xxx,-yyy,'k')
            hold off
            xlabel('x')
            ylabel('r')
            axis equal
            grid on
            axis([-xRange xRange -yRange yRange])
            title(['t=' num2str((i-1)*deltaT)])
            drawnow
                    
            CountSnap = CountSnap+1;
            
        end
            
    end
    
    if plotTracers==1
        
        if plotVelFieldTrace==1
            finess = MeshFine;
            xxx = linspace(-xRange,xRange,finess);
            yyy = linspace(0,yRange,finess);
            [Xfix,Yfix] = meshgrid(xxx,yyy);

            %compute velovity on fixed grid
            [Ufix,Vfix] = velocity_field_extensional_tracers(a,b,risy(:,ite),visc,Ca,Xfix,Yfix,0);
        
        
            %get colorbar limits
            if iteSave==1

                Umax = max(max(Ufix));
                Vmax = max(max(Vfix));
                IntensMax = max(max(sqrt(Ufix.^2+Vfix.^2)));

                Umin = min(min(Ufix));
                Vmin = min(min(Vfix));
                IntensMin = min(min(sqrt(Ufix.^2+Vfix.^2)));

            end
        
        end
    
        %plot tracers
        figure(1)
        %contourf(Xfix,Yfix,sqrt(Ufix.^2+Vfix.^2),10,'LineStyle','none')
        if plotVelFieldTrace==1
            contourf(Xfix,Yfix,Vfix,10,'LineStyle','none')
            hold on
            %contourf(Xfix,-Yfix,sqrt(Ufix.^2+Vfix.^2),10,'LineStyle','none')
            contourf(Xfix,-Yfix,Vfix,10,'LineStyle','none')
        end
        plot(X,-Y,'.b','MarkerSize',10)
        hold on
        plot(X,Y,'.b','MarkerSize',10)
        plot(a,b,'k',a,-b,'k')
        axis equal
        axis([-xRange xRange -yRange yRange])
        grid on
        xlabel('x')
        ylabel('r')
        colorbar
        %caxis([Vmin Vmax])
        hold off
        title('Tracers in extensional flow')
        drawnow

        %compute velovity
        [U,V] = velocity_field_extensional_tracers(a,b,risy(:,ite),visc,Ca,X,Y,0);
        
        %replace particle out
        InOut = FigInOut(a,b,X(:),Y(:));
        InOut = reshape(InOut,numel(y),numel(x));
        X = (InOut>pi).*X + (InOut<pi)*0;
        Y = (InOut>pi).*Y + (InOut<pi)*0;
        U = (InOut>pi).*U + (InOut<pi)*0;
        V = (InOut>pi).*V + (InOut<pi)*0;

        %update tracers
        X = X+deltaT*U*step;
        Y = Y+deltaT*V*step;
        
        if saveTracers==1

            name = 'extensionalTracers';
            print('-dpng','-loose','-r100',[dest name sprintf('%04d',iteSave) '.png'])
            
            
        end
    
    end 
        
    if plotVelField==1
        
        finess = MeshFine;
        x = linspace(-xRange,xRange,finess);
        y = linspace(0,yRange,finess);
        [X,Y] = meshgrid(x,y);
        
        %compute velocity
        [U,V] = velocity_field_extensional_tracers(a,b,risy(:,ite),visc,Ca,X,Y,0);
        
        %get colorbar limits
        if iteSave==1
           
            Umax = max(max(U));
            Vmax = max(max(V));
            IntensMax = max(max(sqrt(U.^2+V.^2)));
            
            Umin = min(min(U));
            Vmin = min(min(V));
            IntensMin = min(min(sqrt(U.^2+V.^2)));
            
        end
        
        %plot vel field
        figure(2)
        contourf(X,Y,V,10,'LineStyle','none')
        hold on
        contourf(X,-Y,V,10,'LineStyle','none')
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

            name = 'extensionalRadialVelField';
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

            name = 'extensionalAxialVelField';
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
        hold on
        contourf(X,-Y,sqrt(U.^2+V.^2),10,'LineStyle','none')
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

            name = 'extensionalIntensityVelField';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',iteSave) '.png'])
            
            
        end
        
    end
    
end