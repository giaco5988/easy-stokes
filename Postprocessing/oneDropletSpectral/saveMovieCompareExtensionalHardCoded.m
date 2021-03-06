%post processing of test time stepping spectral

function saveMovieCompareExtensionalHardCoded(Tstore,Ystore,ElongStore,f1Store,f2Store,color,stepSave,nameNorm,namePhase,nameShape,saveLoop,savePlotEnd,solidLine,PARAM,saveDest,saveMovie)

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',3,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%plot elongation
figure
        
    %get data
    elong2 = ElongStore{2};
    T2 = Tstore{2};
    elong3 = ElongStore{3};
    T3 = Tstore{3};
    elong4 = ElongStore{4};
    T4 = Tstore{4};
    elong5 = ElongStore{5};
    T5 = Tstore{5};
    
    %cut data
    [~,indFinal2] = min(abs(T2-savePlotEnd(2)));
    T2 = T2(1:stepSave(2):indFinal2);
    elong2 = elong2(1:stepSave(2):indFinal2);
    [~,indFinal3] = min(abs(T3-savePlotEnd(3)));
    T3 = T3(1:stepSave(3):indFinal3);
    elong3 = elong3(1:stepSave(3):indFinal3);
    [~,indFinal4] = min(abs(T4-savePlotEnd(4)));
    T4 = T4(1:stepSave(4):indFinal4);
    elong4 = elong4(1:stepSave(4):indFinal4);
    [~,indFinal5] = min(abs(T5-savePlotEnd(5)));
    T5 = T5(1:stepSave(5):indFinal5);
    elong5 = elong5(1:stepSave(5):indFinal5);
        
    count = 0;
    for hhh = 1:max(max(numel(T2,T3),max(numel(T4,T5))))
            
           for k = 1:saveLoop(1)-1
               
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
                axis([0 150 1 10])
                grid on
                hold on
                title('')
                xlabel('t')
                ylabel('L')
               
           end
        
%         if solidLine(i)==1
%             lineType = '-';
%         else
%             lineType = '--';
%         end
            
        display(['Save ' nameNorm ' ' num2str(count) ' of ' num2str(numel(1:max(max(numel(T2,T3),max(numel(T4,T5))))))])
        
        if hhh<=numel(T2)
            %hhh
            semilogy(T2(1:hhh),elong2(1:hhh),'Color',color(2,:))
            hold on
            semilogy(T2(hhh),elong2(hhh),'.','MarkerSize',35,'Color',color(2,:))
        else
            semilogy(T2,elong2,'Color',color(2,:))
            hold on
            semilogy(T2(end),elong2(end),'.','MarkerSize',35,'Color',color(2,:))
        end
        
        if hhh<=numel(T3)
            semilogy(T3(1:hhh),elong3(1:hhh),'Color',color(3,:))
            semilogy(T3(hhh),elong3(hhh),'*','MarkerSize',15,'Color',color(3,:))
        else
            semilogy(T3,elong3,'Color',color(3,:))
            hold on
            semilogy(T3(end),elong3(end),'*','MarkerSize',15,'Color',color(3,:))
        end
        
        if hhh<=numel(T4)
            semilogy(T4(1:hhh),elong4(1:hhh),'Color',color(4,:))
            semilogy(T4(hhh),elong4(hhh),'.','MarkerSize',35,'Color',color(4,:))
        else
            semilogy(T4,elong4,'Color',color(4,:))
            hold on
            semilogy(T4(end),elong4(end),'.','MarkerSize',35,'Color',color(4,:))
        end
        
        if hhh<=numel(T5)
            semilogy(T5(1:hhh),elong5(1:hhh),'Color',color(5,:))
            semilogy(T5(hhh),elong5(hhh),'.','MarkerSize',35,'Color',color(5,:))
        else
            semilogy(T5,elong5,'Color',color(5,:))
            hold on
            semilogy(T5(end),elong5(end),'.','MarkerSize',35,'Color',color(5,:))
        end
        
        axis([0 150 1 10])
        grid on
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

    %plot phase space
    figure
        
    %get data
    f1_2 = f1Store{2};
    f2_2 = f2Store{2};
    f1_3 = f1Store{3};
    f2_3 = f2Store{3};
    f1_4 = f1Store{4};
    f2_4 = f2Store{4};
    f1_5 = f1Store{5};
    f2_5 = f2Store{5};
    
    %cut data
    f1_2 = f1_2(1:stepSave(2):indFinal2);
    f2_2 = f2_2(1:stepSave(2):indFinal2);
    f1_3 = f1_3(1:stepSave(3):indFinal3);
    f2_3 = f2_3(1:stepSave(3):indFinal3);
    f1_4 = f1_4(1:stepSave(4):indFinal4);
    f2_4 = f2_4(1:stepSave(4):indFinal4);
    f1_5 = f1_5(1:stepSave(5):indFinal5);
    f2_5 = f2_5(1:stepSave(5):indFinal5);
        
    count = 0;
    for hhh = 1:1:max(max(numel(f1_2,f1_3),max(numel(f1_4,f1_5))))
            
        
        if hhh<=numel(T2)
            plot(f1_2(1:hhh),f2_2(1:hhh),'Color',color(2,:))
            hold on
            plot(f1_2(hhh),f2_2(hhh),'.','MarkerSize',35,'Color',color(2,:))
        else
            plot(f1_2,f2_2,'Color',color(2,:))
            hold on
            plot(f1_2(end),f2_2(end),'.','MarkerSize',35,'Color',color(2,:))
        end
        
        if hhh<=numel(T3)
            plot(f1_3(1:hhh),f2_3(1:hhh),'Color',color(3,:))
            hold on
            plot(f1_3(hhh),f2_3(hhh),'*','MarkerSize',15,'Color',color(3,:))
        else
            plot(f1_3,f2_3,'Color',color(3,:))
            hold on
            plot(f1_3(end),f2_3(end),'*','MarkerSize',15,'Color',color(3,:))
        end
        
        if hhh<=numel(T4)
            plot(f1_4(1:hhh),f2_4(1:hhh),'Color',color(4,:))
            hold on
            plot(f1_4(hhh),f2_4(hhh),'.','MarkerSize',35,'Color',color(4,:))
        else
            plot(f1_4,f2_4,'Color',color(4,:))
            hold on
            plot(f1_4(end),f2_4(end),'.','MarkerSize',35,'Color',color(4,:))
        end
        
        if hhh<=numel(T5)
            plot(f1_5(1:hhh),f2_5(1:hhh),'Color',color(5,:))
            hold on
            plot(f1_5(hhh),f2_5(hhh),'.','MarkerSize',35,'Color',color(5,:))
        else
            plot(f1_5,f2_5,'Color',color(5,:))
            hold on
            plot(f1_5(end),f2_5(end),'.','MarkerSize',35,'Color',color(5,:))
        end
            
        axis([0 1.2 0 1.2])
        grid on
        display(['Save ' namePhase ' ' num2str(count) ' of ' num2str(numel(1:max(max(numel(f1_2,f1_3),max(numel(f1_4,f1_5))))))])
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

% %plot residuals
% figure
% for i = saveLoop
%         
%     %get data
%     T = Tstore{i};
%     res = resStore{i};
%     
%     %cut data
%     indFinal = manyIndFinal(i);
%     T = T(1:indFinal);
%     res = res(1:indFinal);
%         
%     count = 0;
%     for hhh = 1:stepSave:max(numel(T))
%             
%         %plot previous shapes
%         if i>1
%             
%            for k = 1:i-1
%                
%                 %get data
%                 Tbefore = Tstore{k};
%                 resBefore = resStore{k};
% 
%                 %cut data
%                 indFinal = manyIndFinal(k);
%                 Tbefore = Tbefore(1:indFinal);
%                 resBefore = resBefore(1:indFinal);
%                 
%                 if solidLine(k)==1
%                     lineType = '-';
%                 else
%                     lineType = '--';
%                 end
%                 
%                 semilogy(Tbefore,resBefore,lineType,'Color',color(k,:))
%                 axis([0 250 1e-8 1e2])
%                 %grid on
%                 hold on
%                 title('')
%                 xlabel('f_1')
%                 ylabel('f_2')
%                
%            end
%             
%         end
%         
%         if solidLine(i)==1
%             lineType = '-';
%         else
%             lineType = '--';
%         end
%             
%         display(['Save ' nameRes ' ' num2str(count) ' of ' num2str(numel(1:stepSave:max(numel(T))))])
%         semilogy(T(1:hhh),res(1:hhh),lineType,'Color',color(i,:))
%         axis([0 250 1e-8 1e2])
%         grid on
%         hold on
%         semilogy(T(hhh),res(hhh),'.','MarkerSize',35,'Color',color(i,:))
%         title('')
%         xlabel('t')
%         ylabel('$$ ||(\mathbf{u}-\mathbf{u}_{d}) \cdot \mathbf{n}||_\infty $$','Interpreter','latex')
%         drawnow
%         hold off
%         if saveMovie(3)==1
%             print('-dpng','-loose','-r100',[saveDest nameRes sprintf('%03d',count) '.png'])
%         end
%         
%         count = count+1;
%         
%     end
%     
% end

%figure in pixel
width = 600;
height = 200;

%plot shape
    
    figure
    %fig = gcf;
    %fig.Position = [400 200 width height];
        
    %get data
    Y2 = Ystore{2};
    Y2 = Y2(1:stepSave(2):indFinal2,:);
    Y3 = Ystore{3};
    Y3 = Y3(1:stepSave(3):indFinal3,:);
    Y4 = Ystore{4};
    Y4 = Y4(1:stepSave(4):indFinal4,:);
    Y5 = Ystore{5};
    Y5 = Y5(1:stepSave(5):indFinal5,:);
        
    count = 0;
    for hhh = 1:max(max(numel(T2,T3),max(numel(T4,T5))))
        
        if hhh<=numel(T2)
            xMode2 = Y2(hhh,1:2:end-1);
            yMode2 = Y2(hhh,2:2:end);
            [x2,y2] = fromModesToGrid(xMode2,yMode2,PARAM);
            plot(x2,y2,'Color',color(2,:))
            hold on
            plot(x2,-y2,'Color',color(2,:))
        else
            plot(x2,y2,'Color',color(2,:))
            hold on
            plot(x2,-y2,'Color',color(2,:))
        end
        
        if hhh<=numel(T3)
            xMode3 = Y3(hhh,1:2:end-1);
            yMode3 = Y3(hhh,2:2:end);
            [x3,y3] = fromModesToGrid(xMode3,yMode3,PARAM);
            plot(x3,y3+3,'Color',color(3,:))
            hold on
            plot(x3,-y3+3,'Color',color(3,:))
        else
            plot(x3,y3+3,'Color',color(3,:))
            hold on
            plot(x3,-y3+3,'Color',color(3,:))
        end
        
        if hhh<=numel(T4)
            xMode4 = Y4(hhh,1:2:end-1);
            yMode4 = Y4(hhh,2:2:end);
            [x4,y4] = fromModesToGrid(xMode4,yMode4,PARAM);
            plot(x4,y4+6,'Color',color(4,:))
            hold on
            plot(x4,-y4+6,'Color',color(4,:))
        else
            plot(x4,y4+6,'Color',color(4,:))
            hold on
            plot(x4,-y4+6,'Color',color(4,:))
        end
        
        if hhh<=numel(T5)
            xMode5 = Y5(hhh,1:2:end-1);
            yMode5 = Y5(hhh,2:2:end);
            [x5,y5] = fromModesToGrid(xMode5,yMode5,PARAM);
            plot(x5,y5+9,'Color',color(5,:))
            hold on
            plot(x5,-y5+9,'Color',color(5,:))
        else
            plot(x5,y5+9,'Color',color(5,:))
            hold on
            plot(x5,-y5+9,'Color',color(5,:))
        end
        hold off
        
        %display
        display(['Save ' nameShape ' ' num2str(count) ' of ' num2str(max(max(numel(T2,T3),max(numel(T4,T5)))))])
        axis equal
        axis([-10 10 -1 12])
        grid off
        hold on
        title('')
        axis off
        drawnow
        hold off
        if saveMovie(3)==1
            print('-dpng','-loose','-r100',[saveDest nameShape sprintf('%03d',count) '.png'])
        end
        
        count = count+1;
            
    end
