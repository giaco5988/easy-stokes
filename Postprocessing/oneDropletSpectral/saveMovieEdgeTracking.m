%post processing of test time stepping spectral

function saveMovieEdgeTracking(Tstore,Ystore,ElongStore,f1Store,f2Store,resStore,color,stepSave,nameNorm,namePhase,nameRes,nameShape,saveLoop,savePlotEnd,solidLine,PARAM,saveDest,saveMovie)

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',3,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%plot elongation
figure
manyIndFinal = zeros(numel(saveLoop),1);
for i = saveLoop
        
    %get data
    elong = ElongStore{i};
    T = Tstore{i};
    
    %cut data
    [~,indFinal] = min(abs(T-savePlotEnd(i)));
    T = T(1:indFinal);
    elong = elong(1:indFinal);
    
    manyIndFinal(i) = indFinal;
        
    count = 0;
    for hhh = 1:stepSave:max(numel(T))
            
        %plot previous shapes
        if i>1
            
           for k = 1:i-1
               
                %get data
                elongBefore = ElongStore{k};
                Tbefore = Tstore{k};
                
                %cut data
                [~,indFinal] = min(abs(Tbefore-savePlotEnd(k)));
                Tbefore = Tbefore(1:indFinal);
                elongBefore = elongBefore(1:indFinal);
                
                manyIndFinal(k) = indFinal;
                
                if solidLine(k)==1
                    lineType = '-';
                else
                    lineType = '--';
                end
                
                semilogy(Tbefore,elongBefore,lineType,'Color',color(k,:))
                axis([0 250 1 10])
                grid on
                hold on
                title('')
                xlabel('t')
                ylabel('L')
               
           end
            
        end
        
        if solidLine(i)==1
            lineType = '-';
        else
            lineType = '--';
        end
            
        display(['Save ' nameNorm ' ' num2str(count) ' of ' num2str(numel(1:stepSave:max(numel(T))))])
        semilogy(T(1:hhh),elong(1:hhh),lineType,'Color',color(i,:))
        axis([0 250 1 10])
        grid on
        hold on
        semilogy(T(hhh),elong(hhh),'.','MarkerSize',35,'Color',color(i,:))
        title('')
        xlabel('t')
        ylabel('L')
        drawnow
        hold off
        if saveMovie(1)==1
            print('-dpng','-loose','-r100',[saveDest nameNorm sprintf('%03d',count) '.png'])
        end
        
        count = count+1;
        
    end
    
end

%plot phase space
figure
for i = saveLoop
        
    %get data
    f1 = f1Store{i};
    f2 = f2Store{i};
    
    %cut data
    indFinal = manyIndFinal(i);
    f1 = f1(1:indFinal);
    f2 = f2(1:indFinal);
        
    count = 0;
    for hhh = 1:stepSave:max(numel(f1))
            
        %plot previous shapes
        if i>1
            
           for k = 1:i-1
               
                %get data
                f1Before = f1Store{k};
                f2Before = f2Store{k};

                %cut data
                indFinal = manyIndFinal(k);
                f1Before = f1Before(1:indFinal);
                f2Before = f2Before(1:indFinal);
                
                if solidLine(k)==1
                    lineType = '-';
                else
                    lineType = '--';
                end
                
                plot(f1Before,f2Before,lineType,'Color',color(k,:))
                axis([0.2 1 0 0.8])
                %grid on
                hold on
                title('')
                xlabel('f_2')
                ylabel('f_4')
               
           end
            
        end
        
        if solidLine(i)==1
            lineType = '-';
        else
            lineType = '--';
        end
            
        display(['Save ' namePhase ' ' num2str(count) ' of ' num2str(numel(1:stepSave:max(numel(f1))))])
        plot(f1(1:hhh),f2(1:hhh),lineType,'Color',color(i,:))
        axis([0.2 1 0 0.8])
        grid on
        hold on
        plot(f1(hhh),f2(hhh),'.','MarkerSize',35,'Color',color(i,:))
        title('')
        xlabel('f_2')
        ylabel('f_4')
        drawnow
        hold off
        if saveMovie(2)==1
            print('-dpng','-loose','-r100',[saveDest namePhase sprintf('%03d',count) '.png'])
        end
        
        count = count+1;
        
    end
    
end

%plot residuals
figure
for i = saveLoop
        
    %get data
    T = Tstore{i};
    res = resStore{i};
    
    %cut data
    indFinal = manyIndFinal(i);
    T = T(1:indFinal);
    res = res(1:indFinal);
        
    count = 0;
    for hhh = 1:stepSave:max(numel(T))
            
        %plot previous shapes
        if i>1
            
           for k = 1:i-1
               
                %get data
                Tbefore = Tstore{k};
                resBefore = resStore{k};

                %cut data
                indFinal = manyIndFinal(k);
                Tbefore = Tbefore(1:indFinal);
                resBefore = resBefore(1:indFinal);
                
                if solidLine(k)==1
                    lineType = '-';
                else
                    lineType = '--';
                end
                
                semilogy(Tbefore,resBefore,lineType,'Color',color(k,:))
                axis([0 250 1e-8 1e2])
                %grid on
                hold on
                title('')
                xlabel('f_1')
                ylabel('f_2')
               
           end
            
        end
        
        if solidLine(i)==1
            lineType = '-';
        else
            lineType = '--';
        end
            
        display(['Save ' nameRes ' ' num2str(count) ' of ' num2str(numel(1:stepSave:max(numel(T))))])
        semilogy(T(1:hhh),res(1:hhh),lineType,'Color',color(i,:))
        axis([0 250 1e-8 1e2])
        grid on
        hold on
        semilogy(T(hhh),res(hhh),'.','MarkerSize',35,'Color',color(i,:))
        title('')
        xlabel('t')
        ylabel('$$ ||(\mathbf{u}-\mathbf{u}_{d}) \cdot \mathbf{n}||_\infty $$','Interpreter','latex')
        drawnow
        hold off
        if saveMovie(3)==1
            print('-dpng','-loose','-r100',[saveDest nameRes sprintf('%03d',count) '.png'])
        end
        
        count = count+1;
        
    end
    
end

%figure in pixel
width = 600;
height = 200;

%plot shape
for i = saveLoop
    
    figure
    %fig = gcf;
    %fig.Position = [400 200 width height];
        
    %get data
    Y = Ystore{i};
    indFinal = manyIndFinal(i);
        
    count = 0;
    for hhh = 1:stepSave:indFinal
        
        %plot previous shapes
        if i>1
            
           for k = 1:i-1
               
                shift = (k-1)*3;
               
                %get data
                Ybefore = Ystore{k};
                indFinal = manyIndFinal(k);
                
                %get modes
                xMode = Ybefore(indFinal,1:2:end-1);
                yMode = Ybefore(indFinal,2:2:end);

                %get nodal values
                [x,y] = fromModesToGrid(xMode,yMode,PARAM);
                
                plot(x,y-shift,'Color',color(k,:))
                axis equal
                axis([-10 10 -2 2])
                grid off
                hold on
                plot(x,-y-shift,'Color',color(k,:))
                title('')
                axis off
               
           end
            
        end
        
        shift = (i-1)*3;
        
        %get modes
        xMode = Y(hhh,1:2:end-1);
        yMode = Y(hhh,2:2:end);
        
        %get nodal values
        [x,y] = fromModesToGrid(xMode,yMode,PARAM);
        
        %display
        display(['Save ' nameShape ' ' num2str(count) ' of ' num2str(numel(1:stepSave:indFinal))])
        plot(x,y-shift,'Color',color(i,:))
        axis equal
        axis([-10 10 -12 1])
        grid off
        hold on
        plot(x,-y-shift,'Color',color(i,:))
        title('')
        axis off
        drawnow
        hold off
        if saveMovie(4)==1
            print('-dpng','-loose','-r100',[saveDest nameShape sprintf('%03d',count) '.png'])
        end
        
        count = count+1;
            
    end
    
end
