%plot comparison between linear and non linear simulations at different
%time step

clear variables
close all

%load('~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/flow_field/ellipsoidal/q=200_visc=0.5_Dt=0.001_loop=15000_DELTA=0.05_Ca=6_RK2.mat')

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,...
    'defaultpatchlinewidth',.7);

nFrames = 1000;
step = 10;

scrsz = get(0,'ScreenSize');

%otptions
int = 1;            % plot interface
PlotVelField = 1;   % plot velocity field
CM_vel = 0;         % plot center of mass velocity
save = 0;           % save png
dest = '~/Documents/research_notes/APS/movie/dropOblateProlate/';

%center image
xRange = 3;   yRange = 3;

load('~/Documents/MATLAB/droplet_simulations/results/extensional/viscosity/Extensional_q=100_visc=1_Dt=0.004_loop=10000_Ca=0.08_VR=1_RK2.mat')

%initialization
manyD = zeros(1,nFrames-2);

for i = 3:nFrames
    
    ite = (i-1)*step+1;
    
    disp([num2str(i) ' of ' num2str(nFrames)])
    
    if CM_vel==1
    
        a_before = risa(:,ite-2)';
        b_before = risb(:,ite-2)';
        a = risa(:,ite-1)';
        b = risb(:,ite-1)';

        %figure out center of mass velocity
        xcm_old = center_mass(a_before,b_before);
        xcm = center_mass(a,b);
        v_centermass = -(xcm-xcm_old)/deltaT/checkpoint;

        figure(2)
        hold on
        plot(i*checkpoint*deltaT,v_centermass,'o')
        hold off
        xlabel('t')
        ylabel('vcm')
        grid on
    
    end
    
    if int==1
        
        %load('~/Documents/MATLAB/droplet_simulations/results/rising_dropletOLD/lambda=5_Ca=6/oblate/q=200_visc=5_Dt=0.005_loop=16000_DELTA=-0.14_Ca=6_RK2.mat')

        a = risa(:,ite)';
        b = risb(:,ite)';

        %center of mass
        xcm = center_mass(a,b);

        %compute the spline coeff
        [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);

        %compute splines coordinates
        t = 0:0.05:0.9;
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

        figure(1)
        plot(xxx,yyy,'-k',xxx,-yyy,'-k','LineWidth',2)
        axis equal
        axis off
        hold on
        axis([-xRange xRange -yRange xRange])
        %drawnow
        
        if PlotVelField==1
            
            %compute velocity field
            [U,V,X,Y] = velocity_field_extensional(a,b,risy(:,ite),visc,Ca,xRange,yRange,0);
            
            figure(1)
            contour(X,Y,sqrt(U.^2+V.^2),20)
            hold on
            contour(X,-Y,sqrt(U.^2+V.^2),20)
            axis equal
            
            %plot velocity field
            quiver(X,Y,U,V,'r')
            quiver(X,-Y,U,-V,'r')

            plot(a,b,'k-',a,-b,'k-','LineWidth',2)
            plot(a,b,'k-',a,-b,'k-','LineWidth',2)
            drawnow
            
        end
        
        %compute ellipticity
        L = max(xxx)-min(xxx);  B = 2*max(yyy);
        D = (L-B)/(L+B);
        manyD(i-2) = D;

        if save==1

            name = 'extensional';
            print('-dpng','-loose','-r100',[dest name sprintf('%03d',i-2) '.png'])
    %         set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
    %         saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/interface' num2str(i) '.png'])
        end

        hold off

    end
    
    
end

figure
plot(3:nFrames,manyD)
xlabel('iterations')
ylabel('D')
grid on



