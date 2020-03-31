%post processing of bifurcation in extensional flow

clear variables
close all

%parameters for upload
% manyN = [100 100 100 100 50 100 50];
% manyLambda = [0.02 0.05 0.1 0.5 1 5 10];
% manyCaUp = [0.4 0.4 0.4 0.4 0.2 0.4 0.2];
manyN = [50 50 50 50 50 50];
manyLambda = [0.1 0.3 1 10];
manyCaUp = [0 0 0 0 0 0 0 0 0 0 0];
CaDown = -1;
BC = 1;
legendre = 1;
resUp = 0;

%options
plotMaxCa = 0;  plotMinCa = 1;
plotMaxD = 0;
plotCutD = 0;   CutD = 2;
getIte = 1;
plotTipCurvature = 0;
checkVolume = 0;
plotCurv = 0;
rescaleShape = 0;   rescaleLaw = 2;
plotEigWithColors = 0;  line = 1;
chooseBifParameter = 2;  % 1 is the classical D, 2 id the max(x), 3 is max(x) - position of the center of mass, 4 if the max(abs(x)), 5 is surface area
plotCacritHalf = 0;

%initialize
maxDeformationPar = zeros(numel(manyLambda),1);
maxCa = zeros(numel(manyLambda),1);
minCa = zeros(numel(manyLambda),1);
DmaxCa = zeros(numel(manyLambda),1);
Dcut = zeros(numel(manyLambda),1);
cellLegend = cell(size(manyLambda));
shiftHalf = 0;
Lup = zeros(numel(manyLambda),1);
Ldown = zeros(numel(manyLambda),1);

for k = 1:numel(manyLambda)
    
    display([num2str(k) ' of ' num2str(numel(manyLambda))])
    
    %current parameter
    lambda = manyLambda(k);
    n = manyN(k);
    CaUp = manyCaUp(k);
    cellLegend(k) = {['\lambda=' num2str(lambda)]};

    %load data
    if resUp==1
        res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/';
    elseif resUp==0
        res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/server/';
    elseif resUp==2
        res = '/Users/Giacomo/Documents/MATLAB/droplet_simulations/results/edgeTracking/forPaper/data/';
    end
    filename = [res 'newtonMethodSpectralXYmodesCont_n=' num2str(n) '_CaUp=' num2str(CaUp) '_CaDown='  num2str(CaDown) '_visc=' num2str(lambda) '_Legendre=' num2str(legendre) '_BC=' num2str(BC) '.mat'];
    load(filename)

    %for which point I plot
    Dplot = [];
    %Dplot = [2.7 3.01 3.026 3.136 3.14];
    %Dplot = [2.7 3.01 3.026 3.136 3.14];
    CaPlot = [0.1015 0.093 0.097 0.095 0.0965 0.09555];

    %CaPlot = [0.05 0.07 0.1];
    CaPlot = 0.03;
    %Dplot = [1.07 1.13 1.2 1.36];
    %Dplot = [2 2.26 2.41 2.19];
    %Dplot = 1.07;
    %Dplot = 2;

    ind0 = find(manyD==0,1)-getIte;
    manyD = manyD(1:ind0);
    manyCa = manyCa(1:ind0);
    %unstEigen = 1;
    unstEigen = unstEigen(1:ind0)+1;
    manyA = manyA(:,1:ind0);
    manyB = manyB(:,1:ind0);
    errV = zeros(ind0,1);
    
    %try scaling
    [maxDeformationPar(k),ind] = max(manyD);
    [~,cut] = min(abs(manyD(1:ind)-CutD));
    Dcut(k) = manyCa(cut);
    %maxDormationParCa(k) = manyCa(inD);
    [maxCa(k),indMaxCa] = max(manyCa);
    minCa(k) = min(manyCa);

    %compute bifurcation parameter
    if chooseBifParameter==1

        L = manyA(1,:)-manyA(end,:);
        B = 2*manyB(round(PARAM.n/2)+1,:);
        %B = 2*max(manyB);
        manyD = (L-B)./(L+B);
        figure(1)
        ylabel('D')

    elseif chooseBifParameter==2

        manyD = manyA(1,:);
        figure(1);
        ylabel('max(x)')

    elseif chooseBifParameter==3

        for i = 1:ind0

            xxx = manyA(:,i);   yyy = manyB(:,i);
            manyD(i) = max(manyA(:,i)) - CenterMassCurvAxisSpectral(xxx,yyy,PARAM);
            %manyD(i) = CenterMassCurvAxisSpectral(xxx,yyy,PARAM);

        end

        figure(1)
        %subplot(2,2,1)
        ylabel('||x||_{\infty}-x_{cm}')

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

    end
    DmaxCa(k) = manyD(indMaxCa);

    %bifurcation diagram
    figure(1);
    hold on
    plot(manyCa,manyD)
    %plot(manyCa,manyD)
    xlabel('Ca')
    title('Bifurcation diagram')
    grid on

    col = 'krgmby';
    col = repmat(col,1,10);
    
    if plotCacritHalf == 1 && k>2
       
        [~,indDown] = min(abs(manyCa(1:indMaxCa)-maxCa(k)*0.5));
        [~,indUp] = min(abs(manyCa(indMaxCa:end)-maxCa(k)*0.5));
        %[~,indDown] = min(abs(manyCa(1:indMaxCa)-0.05));
        %[~,indUp] = min(abs(manyCa(indMaxCa:end)-0.05));
        xDown = manyA(:,indDown);
        yDown = manyB(:,indDown);
        xUp = manyA(:,indUp+indMaxCa);
        yUp = manyB(:,indUp+indMaxCa);
        plot(0.5*[maxCa(k) maxCa(k)],[manyD(:,indDown) manyD(:,indUp+indMaxCa-1);],'o')
        
        figure(18)
        hold on
        plot(xDown+shiftHalf,yDown,'k')
        plot(xDown+shiftHalf,-yDown,'k')
        plot(xUp+shiftHalf,yUp+3,'k')
        plot(xUp+shiftHalf,-yUp+3,'k')
        axis off
        axis equal
        
        shiftHalf = shiftHalf+8;
        Lup(k) = manyD(:,indUp+indMaxCa-1);
        Ldown(k) = manyD(:,indDown);
        
    elseif plotCacritHalf == 1 && k<=2
        
        maxCa(k) = 0.2;
        [~,indDown] = min(abs(manyCa-maxCa(k)));
        xDown = manyA(:,indDown);
        yDown = manyB(:,indDown);
        plot(maxCa(k),manyD(:,indDown),'o')
        
        figure(18)
        hold on
        plot(xDown+shiftHalf,yDown,'k')
        plot(xDown+shiftHalf,-yDown,'k')
        axis off
        axis equal
        
        shiftHalf = shiftHalf+8;
        Ldown(k) = manyD(:,indDown);
        
    end

    if plotEigWithColors==1

        for i = min(unstEigen(1:ind0)):numel(col)

            if i>min(unstEigen(1:ind0)); hold on; end

            indCol = find(unstEigen==i,numel(unstEigen),'first');

            figure(1)
            hold on
            %subplot(2,2,1)
            if line==0
                colHere = [col(i) 'o'];
            elseif line==1
                colHere = [col(i) '-'];
                hold on
            end
            plot(manyCa(indCol),manyD(indCol),colHere)
            %hold off
            title(['Bifurcation diagram \lambda=' num2str(lambda)])
            grid on

        end

    end

    if numel(Dplot)>=1

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
            shift = 2*countShapes;
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
            xlabel('x')
            ylabel('r')
            title('Steady states')

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
        plot(PARAM.t,K1,'-')
        hold on
        plot(PARAM.t,K2,'-')
        grid on
        xlabel('t')
        ylabel('K')
        legend('K_1','K_2','Location','Best')
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

end

figure(1)
legend(cellLegend,'Location','Best')
plot(maxCa,DmaxCa,'ko')

if plotMaxD==1

    figure
    loglog(manyLambda,maxDeformationPar/maxDeformationPar(1),'o-')
    grid on
    xlabel('\lambda')
    ylabel('max(D)')

end

if plotMaxCa==1

    figure
    loglog(manyLambda,maxCa,'o-')
    hold on
    %loglog(manyLambda,manyLambda.^(-1/6)*0.1,'k-')
    ax = loglog(manyLambda,manyLambda.^(-1/6)*0.148,'k-');
    %axis([1e-2 1 0.12 max(manyLambda.^(-1/6)*0.148)])
    axis([1e-2 1 0.12 0.28])
    grid on
    xlabel('\lambda')
    ylabel('max(Ca)')
    legend('Data from simulation','0.148\lambda^{-1/6}','Location','Best')
    set(ax,'YTick',0.12:0.04:0.28)

end

if plotMinCa==1

    figure
    loglog(manyLambda,minCa,'o-')
    hold on
    %loglog(manyLambda,manyLambda.^(-1/6)*0.1,'k-')
    %ax = loglog(manyLambda,manyLambda.^(-1/6)*0.148,'k-');
    %axis([1e-2 1 0.12 max(manyLambda.^(-1/6)*0.148)])
    %axis([1e-2 1 0.12 0.28])
    grid on
    xlabel('\lambda')
    ylabel('min(Ca)')
    %legend('Data from simulation','0.148\lambda^{-1/6}','Location','Best')
    %set(ax,'YTick',0.12:0.04:0.28)

end

if plotCutD==1

    figure
    loglog(manyLambda,Dcut,'o-')
    grid on
    xlabel('\lambda')
    ylabel('max(Ca)')
    title(['Ca for D=' num2str(CutD)])

end














