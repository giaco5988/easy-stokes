%plot droplet curvature

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

here = pwd;
path = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/relaxation/';

%simulation parameters
D = 0.1*(2:7);
%D = 0.8;
visc = 1;
elem = 30;
dt = 0.01;
loop = 1500;
modes = 4;
cellLegend = cell(size(D));

%option
plotShape = 0;

%initialize
kUp = zeros(loop,1);
kSide = zeros(loop,1);
k1Up = zeros(loop,1);
k1Side = zeros(loop,1);
f1 = zeros(loop,1);
f2 = zeros(loop,1);
f3 = zeros(loop,1);

%loop on ellipcisity parameter
for k = 1:numel(D)
    
    %store element in cell for legend
    cellLegend(k) = {num2str(D(k))};
    
    %loop number
    display([num2str(k) ' of ' num2str(numel(D))])

    %loop on simulation step
    for i = 1:loop
        
        display([num2str(i) ' of ' num2str(loop)])
        
        %load data
        filename = [path 'Relaxation_q=' num2str(elem) '_visc=' num2str(visc) '_Dt=' num2str(dt) '_loop=' num2str(loop) '_DELTA=' num2str(D(k)) '_Ca=1_RK2.mat'];
        load(filename)
        
        %current coordinates
        a = risa(:,i)'; b = risb(:,i)';
        
        V = axis_int_gauss_vect(a,b);
        
        %compute spline cofficients
        [ax, bx, cx, dx, ay, by, cy, dy] = spline_symmetric (a, b);
        k1 = curv_spline2(bx,by,cx,cy,dx,dy);

        %save curvature
        kUp(i) = KKK(elem/2+1,i);
        kSide(i) = KKK(1,i);
        k1Up(i) = k1(elem/2+1);
        k1Side(i) = k1(1);

        %plot shape
        if plotShape==1
            figure(1)
            plot([a flip(a)],[b -flip(b)])
            axis equal
            axis([-2 2 -1.2 1.2])
            grid on
            xlabel('z')
            ylabel('r')
            drawnow
        end
        
        %compute legendre modes
        [f,err] = LegendreSerie(a',b',modes,1,V);
        f1(i) = f(2);
        f2(i) = f(3);
        f3(i) = f(4);

    end
    
    figure(2)
    hold on
    plot(kUp,kUp+kSide)
    xlabel('k_s')
    ylabel('k_s+k_e')
    hold off
    grid on
    drawnow
    
    figure(3)
    hold on
    plot(k1Up,k1Up+k1Side)
    xlabel('k1_s')
    ylabel('k1_s+k1_e')
    hold off
    grid on
    drawnow
    
    figure(4)
    hold on
    plot(f1,f2,'-')
    xlabel('f_1')
    ylabel('f_2')
    hold off
    grid on
    drawnow
    
    figure(5)
    if k>1
        hold on
    end
    plot3(f1,f2,f3,'-')
    xlabel('f_1')
    ylabel('f_2')
    zlabel('f_3')
    hold off
    grid on
    drawnow

end


%plot ellispoidal relaxation curve
figure(3)
%compute analytically
D = 0.02:0.01:0.8;
b = ((1.0-D)./(1.0+D)).^(1/3);    %minor axis
a = 1.0./b.^2;                    %major axis
KeKs = b./a.^2 .* (1+a.^3./b.^3);

%compute modes for ellipses
theta = linspace(0,pi,50);
f1Ellipse = zeros(1,numel(D));
f2Ellipse = zeros(1,numel(D));
f3Ellipse = zeros(1,numel(D));
for i = 1:numel(D)
    
    %ellpse volume
    V0 = 4/3*pi;
    
    %major and minor axis
    aHere = a(i);
    bHere = b(i);
    
    %ellipse shape
    x = aHere*cos(theta);
    y = bHere*sin(theta);
    
    %compute legendre
    [fEllipse,err] = LegendreSerie(x',y',modes,1,V0);
    f1Ellipse(i) = fEllipse(2);
    f2Ellipse(i) = fEllipse(3);
    f3Ellipse(i) = fEllipse(4);
    
end

figure(3)
hold on
plot(b./a.^2,KeKs,'k')    %analytical
hold off
grid on
xlabel('K_s')
ylabel('K_s+K_e')
title('Meridional plane')
legend(cellLegend,'Location','Best')

figure(2)
hold on
KeKs = (b./a.^2 + 1./b) .* (1 + 2*a./b.^2./(b./a.^2 + 1./b));
plot(b./a.^2 + 1./b,KeKs,'k')    %analytical
hold off
grid on
xlabel('K_s')
ylabel('K_s+K_e')
title('Full curvature')
legend(cellLegend,'Location','Best')

figure(4)
hold on
plot(f1Ellipse,f2Ellipse,'k')    %analytical
hold off
grid on
xlabel('f_1')
ylabel('f_2')
title('Legendre modes')
legend(cellLegend,'Location','Best')

figure(5)
hold on
plot3(f1Ellipse,f2Ellipse,f3Ellipse,'k')    %analytical
hold off
grid on
xlabel('f_1')
ylabel('f_2')
zlabel('f_3')
title('Legendre modes')
legend(cellLegend,'Location','Best')





