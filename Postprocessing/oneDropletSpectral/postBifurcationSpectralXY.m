%post processing of bifurcation in extensional flow

clear variables
%close all

%parameters for upload
n = 50;
lambda = 1;
CaUp = 1;
CaDown = -1;
BC = 1;
legendre = 1;
res = 1;

%options
subPLOT = 0;
getIte = 230;
plotTipCurvature = 0;
checkVolume = 1;
plotCurv = 0;
rescaleShape = 0;   rescaleLaw = 2;
plotEigWithColors = 0;  line = 1;
chooseBifParameter = 2;  % 1 is the classical D, 2 is the half-width, 3 is max(x) - position of the center of mass, 4 if the max(abs(x)), 5 is surface area, 6 is the center of mass position
plotManyCurv = 1;

%load data
if res==1
    res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/';
elseif res==0
    res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/server/';
end
filename = [res 'newtonMethodSpectralXYmodesCont_n=' num2str(n) '_CaUp=' num2str(CaUp) '_CaDown='  num2str(CaDown) '_visc=' num2str(lambda) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '.mat'];
load(filename)

%for which point I plot
Dplot = [];
Dplot = 2.2;
%Dplot = [0.68 0.35 0.088];

CaPlot = 0.09;
%CaPlot = [-0.2 -0.38 -0.38];

ind0 = find(manyD==0,1)-getIte;
manyD = manyD(1:ind0);
manyCa = manyCa(1:ind0);
%unstEigen = 1;
unstEigen = unstEigen(1:ind0)+1;
manyA = manyA(:,1:ind0);
manyB = manyB(:,1:ind0);
errV = zeros(ind0,1);
tipK1 = zeros(ind0,1);
tipK2 = zeros(ind0,1);
tipK = zeros(ind0,1);
midK1 = zeros(ind0,1);
midK2 = zeros(ind0,1);
midK = zeros(ind0,1);

%compute bifurcation parameter
if chooseBifParameter==1
   
    L = manyA(1,:)-manyA(end,:);
    B = 2*manyB(round(PARAM.n/2)+1,:);
    %B = 2*max(manyB);
    manyD = (L-B)./(L+B);
    figure(1)
    %ylabel('D')
    xyLabelTex('\rm{Ca}','D')
    
elseif chooseBifParameter==2
    
    manyD = (manyA(1,:)-manyA(end,:))/2;
    %manyD = manyA(1,:);
    
    if subPLOT==0
        figure(1);
    elseif subPLOT==1
        subplot(2,2,1)
    end
    %ylabel('L')
    xyLabelTex('\rm{Ca}','L')
    
elseif chooseBifParameter==3
    
    for i = 1:ind0
        
        xxx = manyA(:,i);   yyy = manyB(:,i);
        manyD(i) = max(manyA(:,i)) - CenterMassCurvAxisSpectral(xxx,yyy,PARAM);
        %manyD(i) = CenterMassCurvAxisSpectral(xxx,yyy,PARAM);
        
    end
    
    figure(1)
    %subplot(2,2,1)
    ylabel('||x||_{\infty}-x_{cm}')
    %xyLabelTex('\rm{Ca}','D')
    
elseif chooseBifParameter==4
    
    manyD = max(abs(manyA));
    
    figure(1)
    %subplot(2,2,1)
    ylabel('||x||_{\infty}')
    
elseif chooseBifParameter==5
    
    for i = 1:ind0
        
        xxx = manyA(:,i);   yyy = manyB(:,i);
        manyD(i) = surfaceCurvilinearAxisSpectral(xxx,yyy,PARAM) - 4*pi;
        
    end
    
    figure(1)
    %subplot(2,2,1)
    ylabel('||\Delta S||')
    xyLabelTex('\rm{Ca}','\Delta S')
    
elseif chooseBifParameter==6
    
    for i = 1:ind0
        
        xxx = manyA(:,i);   yyy = manyB(:,i);
        manyD(i) = CenterMassCurvAxisSpectral(xxx,yyy,PARAM);
        
    end
    
    figure(1)
    %subplot(2,2,1)
    %ylabel('z_{\rm{cm}}')
    xyLabelTex('\rm{Ca}','z_{\rm{cm}}')
    
end

%bifurcation diagram
if subPLOT==0
        figure(1);
elseif subPLOT==1
        subplot(2,2,1)
end
hold on
plot(manyCa,manyD,'k')
%plot(manyCa,manyD)
xlabel('Ca')
title(['Bifurcation diagram \lambda=' num2str(lambda)])
grid on

col = 'krgmby';
col = repmat(col,1,10);

if subPLOT==0
    figure
end
if plotEigWithColors==1
    
    for i = min(unstEigen(1:ind0)):numel(col)
        
        if i>min(unstEigen(1:ind0)); hold on; end
        
        indCol = find(unstEigen==i,numel(unstEigen),'first');

        figure(1)
        %subplot(2,2,1)
        if line==0
            colHere = [col(i) 'o'];
        elseif line==1
            colHere = [col(i) '-'];
            hold on
        end
        plot(manyCa(indCol),manyD(indCol),colHere)
        hold off
        title(['Bifurcation diagram \lambda=' num2str(lambda)])
        grid on
    
    end
    
end

if isempty(Dplot)
    
    %plot last shape
    figure(2)
    %subplot(2,2,3)
    x = manyA(:,ind0);   y = manyB(:,ind0);
    plot(x,y,'k')
    hold on
    plot(x,-y,'k')
    axis equal
    grid on
    xlabel('x')
    ylabel('r')
    %title(['Shape at iteration ' num2str(ind0)])
    title(['Shape at Ca=' num2str(manyCa(ind0)) ' D=' num2str(manyD(ind0))])
    axis([-4 4 -2 2])
    
    %subplot(2,2,1)
    figure(1)
    hold on
    plot(manyCa(ind0),manyD(ind0),'.k','MarkerSize',40)
    hold off
    
else

    countShapes = 0;
    K1start = zeros(numel(Dplot),1);
    K1end = zeros(numel(Dplot),1);
    xMinCurv = zeros(numel(Dplot),1);
    yMinCurv = zeros(numel(Dplot),1);
    for i = 1:numel(Dplot)

        display([num2str(i) ' of ' num2str(numel(Dplot))])
        
        if numel(Dplot)==1
            col = 'k';
        else
            col = '';
        end

        [~,ind] = min(sqrt(abs(manyD'-Dplot(i)).^2+abs(manyCa-CaPlot(i)).^2));
        x = manyA(:,ind);   y = manyB(:,ind);

        figure(1)
        %subplot(2,2,1)
        hold on
        plot(manyCa(ind),manyD(ind),['.' col],'MarkerSize',40)
        hold off
        
        %shift sucessive shape
        shift = 4*countShapes;
        %shift = 0;

        figure(2)
        %subplot(3,2,i)
        hold on
        plot([x; flip(x)],[y; -flip(y)]-shift,['-' col])
        grid on
        %axis([max(x)-0.2 max(x)+0.2 -0.2 0.2])
        axis equal
        axis off
        hold off
        %xlabel('x')
        %ylabel('r')
        %title('Steady states')
        
        if plotTipCurvature==1
            
            %compute rho in symmetry axis
            D1 = PARAM.D1;
            D2 = PARAM.D2;

            %compute geomtrical derivaties
            xp = D1*x;    yp = D1*y;
            xpp = D2*x;    ypp = D2*y;

            %compute normal vector
            h = (xp.^2+yp.^2).^(0.5);
            %nx = yp./h;
            ny = -xp./h;

            %compute curvature in meridional plane
            K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);

            %compute curvature in aximuthal direction
            K2 = ny./y;
            
            K1start(i) = K1(1);
            K1end(i) = K1(end);
            
            %take the first convexity starting from the left
            [~,indCurv] = min(K1);
            xMinCurv(i) = max(x)-x(indCurv);
            yMinCurv(i) = y(indCurv);
            
        end
        
        countShapes = countShapes+1;

    end
    
    if plotTipCurvature==1
        
        figure(10)
        semilogy(K1start,'o')
        %semilogy(Dplot,K1start)
        grid on
        xlabel('shape')
        ylabel('K1')
        title('Tips curvature')
       
        hold on
        semilogy(K1end,'o')
        %semilogy(Dplot,K1end)
        grid on
        xlabel('shape')
        ylabel('K1_{last}')
        legend('Left tip','Right tip','Location','Best')
        
        figure(11)
        semilogy(xMinCurv,'o')
        hold on
        semilogy(yMinCurv,'o')
        grid on
        xlabel('shape')
        ylabel('K1')
        title('First cap coordinates')
        legend('x coord','y coord','Location','Best')
        
    end
    
    %rescale and overlap shapes
    if rescaleShape==1
        
        %define rescale
        fitExp = fit((1:numel(Dplot))',K1start,'exp1');
        
        if rescaleLaw==1
            %rescale with curvature
            xStar = K1start./K1start(1);
            yStar = K1start./K1start(1);
        elseif rescaleLaw==2
            %rescale first neck
            xStar = xMinCurv(1)./xMinCurv;
            yStar = yMinCurv(1)./yMinCurv;
        else
            %do nor rescale
            xStar = ones(numel(Dplot),1);
            yStar = ones(numel(Dplot),1);
        end
        
        fitExp2 = fit((1:numel(Dplot))',xMinCurv,'exp1');
        fitExp3 = fit((1:numel(Dplot))',yMinCurv,'exp1');
        
        figure(10)
        semilogy(fitExp.a*exp(fitExp.b*(1:numel(Dplot))),'k')
        legend('Left tip','Right tip','exponential fitting','Location','Best')
        
        figure(11)
        semilogy(fitExp2.a*exp(fitExp2.b*(1:numel(Dplot))),'k')
        semilogy(fitExp3.a*exp(fitExp3.b*(1:numel(Dplot))),'k')
        legend('x coord','y coord','exp fit x','exp fit y','Location','Best')
        
        for i = 1:numel(Dplot)
            
            if numel(Dplot)==1
                col = 'k';
            else
                col = '';
            end

            [~,ind] = min(sqrt(abs(manyD'-Dplot(i)).^2+abs(manyCa-CaPlot(i)).^2));
            x = manyA(:,ind);   y = manyB(:,ind);
            
            if i==1
                
                x = x*xStar(i);
                y = y*yStar(i);
                
                xRef = max(x);
                figure(3)
                
            else
                
                x = x*xStar(i);
                y = y*yStar(i);
                x = x-max(x)+xRef;
                figure(3)
                hold on
                
            end
            
            plot([x; flip(x)],[y; -flip(y)],['-' col])
            grid on
            axis equal
            hold off
            xlabel('x')
            ylabel('r')
            %axis([1.8 2.8 -0.4 0.4])
            title('Rescaled shape')
        
            %width
            width(i) = x(1)-x(end);
            [~,ind] = min(abs(x-2.45));
            yLaw(i) = y(ind);
            
        end
        
        figure(4)
        semilogy(1:numel(Dplot),width,'-o')
        xlabel('shape')
        ylabel('\Delta x')
        title('Droplet width')
        grid on
        
        figure(5)
        plot(1:numel(Dplot),yLaw,'-o')
        xlabel('shape')
        ylabel('y')
        title('y scale')
        grid on
            
    end

end



if plotCurv==1
    
    %coordinates
    x = manyA(:,ind0);   y = manyB(:,ind0);
    
    %compute rho in symmetry axis
    D1 = PARAM.D1;
    D2 = PARAM.D2;
        
    %compute geomtrical derivaties
    xp = D1*x;    yp = D1*y;
    xpp = D2*x;    ypp = D2*y;

    %compute normal vector
    h = (xp.^2+yp.^2).^(0.5);
    %nx = yp./h;
    ny = -xp./h;
        
    %compute curvature in meridional plane
    K1 = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(1.5);
    
    %compute curvature in aximuthal direction
    K2 = ny./y;
    
    %plot last r(theta)
    figure
    %subplot(2,2,4)
    plot(x,K1+K2,'-k')
    %hold on
    %plot(PARAM.t,K2,'-')
    grid on
    xlabel('t')
    ylabel('K')
    %legend('K_1','K_2','Location','Best')
    xyLabelTex('z','K')
    title(['Curvature at Ca=' num2str(manyCa(ind0)) ' D=' num2str(manyD(ind0))])
    
end

%chechVolume
if checkVolume==1
    
   for i = 1:ind0
      
       xxx = manyA(:,i);   yyy = manyB(:,i);
       errV(i) = abs(VolumeCurvilinearAxisSpectral(xxx,yyy,PARAM)-V0)/V0;
       
   end
   
   %subplot(2,2,2)
   figure
   title('error on volume')
   semilogy(errV,'-o')
   xlabel('ite')
   ylabel('errV')
   grid on
    
end

if plotManyCurv==1
    
    for i = 1:ind0
        x = manyA(:,i);
        y = manyB(:,i);
        [K1,K2,K] = curvatureAxisSpectralXY(x,y,PARAM);
        %[K1,K2] = fromModesToGrid(K1,K2,PARAM);
        %K = K1+K2;
        %K = interpLegendreZeroOne([0 0.5],K);
        nnn = numel(x);
        tipK1(i) = K1(1);
        tipK2(i) = K2(1);
        tipK(i) = K(1);
        midK1(i) = K1(nnn/2+1);
        midK2(i) = K2(nnn/2+1);
        midK(i) = K(nnn/2+1);
    end
    figure
    plot(manyCa,tipK1)
    hold on
    plot(manyCa,midK1)
    xlabel('Ca')
    ylabel('K_1')
    legend('Tip curvature','Middle curvature','Location','Best')
    grid on
    
    figure
    plot(manyCa,tipK2)
    hold on
    plot(manyCa,midK2)
    xlabel('Ca')
    ylabel('K_2')
    legend('Tip curvature','Middle curvature','Location','Best')
    grid on
    
    figure
    plot(manyCa,tipK)
    hold on
    plot(manyCa,midK)
    xlabel('Ca')
    ylabel('K')
    legend('Tip curvature','Middle curvature','Location','Best')
    grid on
end

