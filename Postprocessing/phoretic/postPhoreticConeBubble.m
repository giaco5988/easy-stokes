%post processing of a cone moving because of phoretic effects

clear variables
close all

%parameters
manyN = 50*ones(1,1);
nBubble = 1;
%manyN = 1126*ones(1,1);
%manyN = 2141*ones(1,1);
%manyXi = linspace(2,20,10);
manyXi = 10;
manyH = 0.2*ones(1,101);
%manyTheta = 0;
%manyTheta = linspace(0,pi/4,101);
manyTheta = 2/180*pi;
%flux = [0 0 1 0];
M = 1;
res = 1;
rBubble = 1.1;
InPos = 9.9;
BC = 15;

%options
rangeX = 7; rangeY = 3;  shift = 5;
plotSnaphot = 1;  thetaPlot = manyTheta;
threeDup = 0; theta3Dup = {[0 pi] [0 2*pi]};    color3D = {10 [1 1 1]};   transp = [1 0.5];
plotConc = 1;   plotConcLine = 1;   substituiteConc = 1;    coeffSub = 1;
plotCriticalRadius = 1;
plotSutfaceSol = 1;
plotVelField = 0;   labFrame = 1;   % 0 is in motor frame, 1 is in lab frame
checkDecay = 0;     xDecay = logspace(1,5,100); yDecay = zeros(1,numel(xDecay));

%grid for field quantities
nGrid =  10;
x = linspace(-rangeX+shift,rangeX+shift,2*rangeX*nGrid);
y = linspace(0,rangeY,rangeY*nGrid);
[X,Y] = meshgrid(x,y);

%upload directory
if res==1
   dir = '~/Documents/MATLAB/phoreticSwimmer/results/'; 
end

%initialize
U0 = zeros(numel(manyN),1);
U0plot = zeros(numel(manyN),1);
cellLegend = cell(numel(manyN),1);

%countPlot = 1;
doNotPlot = 0;
maxC = [];
for kkk = 1:numel(manyTheta)
    theta = manyTheta(kkk);
    countPlot = 1;
for iii = 1:numel(manyN)
    
    display([num2str(iii*kkk) ' of ' num2str(numel(manyN)*numel(manyTheta))])
    
   %parameters
   h = manyH(iii);
   L = manyXi(iii);
   n = manyN(iii);
   
   %filename
   filename = ['coneWithBubble_n=' num2str(n) '_nBubble=' num2str(nBubble) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '_rBubble=' num2str(rBubble) '_InPos=' num2str(InPos) '.mat'];
   
   %upload data
   %pathUp = [dir ['coneWithBubble_n=' num2str(n) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '_rBubble=' num2str(rBubble) '_InPos=' num2str(InPos) '.mat']];
   pathUp = [dir filename];
   load(pathUp)
   
   %U0(i) = yStokes(end);
            if checkDecay==1 && numel(manyN)==1 && numel(manyTheta)==1

                figure(1)
                if countPlot>1
                    hold on
                end
                [~,~,uField,vField] = computeVelPressField(xDecay,yDecay,x,y,yStokes(1:end-1),0,PARAM);
                Uabs = sqrt(uField.^2+vField.^2);
                loglog(xDecay,Uabs,'o-')
                xlabel('d')
                ylabel('|U|')
                title('Front velocity')
                grid on
                cellLegend{countPlot} = ['\theta=' num2str(manyTheta(iii))];

            end
            
           if plotConc==1 && numel(manyN)==1 && numel(manyTheta)==1
                
               %compute concentration
               figure
               [X,Y,PHIfield] = computeConcentrationField(X,Y,x,y,conc,PARAM,substituiteConc,coeffSub);
               [~,h1] = contourf(X,Y,PHIfield,500);
               hold on
               [~,h2] = contourf(X,-Y,PHIfield,500);
               grid off
               colorbar
               title('Concentration field')
               
               set(h1,'LineColor','none')
               set(h2,'LineColor','none')

                
           end

           frame = '';
           if plotVelField==1 && numel(manyN)==1 && numel(manyTheta)==1
                
               %compute concentration
               figure
               [X,Y,uField,vField] = computeVelPressField(X,Y,x,y,yStokes(1:end-1),0,PARAM);
               if labFrame==0
                    uField = uField-yStokes(end);
                    frame = 'Motor frame';
               elseif labFrame==1
                    frame = 'Lab frame';
               end
               %contourf(X,Y,sqrt(uField.^2+vField.^2))
               %hold on
               %contourf(X,-Y,sqrt(uField.^2+vField.^2))
               %colorbar
               streamslice(X,Y,uField,vField)
               hold on
               streamslice(X,-Y,uField,-vField)
               drawnow
               xlabel('x')
               ylabel('r')
               title('Velocity field')
                
           end
           
           %plot cone
           if numel(manyN)==1 && numel(manyTheta)==1
               plotGeometry(x,y,threeDup,theta3Dup,color3D,transp,PARAM);
               grid on
               axis equal
               axis([-rangeX+shift rangeX+shift -rangeY rangeY])
               xlabel('z')
               ylabel('r')
               axis off
           end
           %title([frame ' \theta=' num2str(theta) ' U=' num2str(yStokes(end))])

           if plotSutfaceSol==1 && numel(manyN)==1 && numel(manyTheta)==1

                %figure
                lRough = [0; cumsum(sqrt(diff(Xsing).^2+diff(Ysing).^2))];
%                 plot(lRough(1:sum(PARAM.n(1:4))),conc(1:sum(PARAM.n(1:4))))
%                 hold on
%                 plot(lRough(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),vSlip)
%                 xlabel('l')
%                 ylabel('c,\nabla c')
%                 grid on
                %title(['\theta=' num2str(theta) ' U=' num2str(yStokes(end))])
                
                figure
                if numel(PARAM.n)==5
                    plot(lRough(sum(PARAM.n(1:4))+1:end)-lRough(sum(PARAM.n(1:4))+1),conc(sum(PARAM.n(1:4))+1:end))
                    title('Flux to bubble')
                elseif numel(PARAM.n)==6
                    subplot(1,2,1)
                    plot(lRough(sum(PARAM.n(1:4))+1:(sum(PARAM.n(1:5))))-lRough(sum(PARAM.n(1:4))+1),conc(sum(PARAM.n(1:4))+1:(sum(PARAM.n(1:5)))))
                    xlabel('l')
                    %legend('flux to bubble 1','Location','Best')
                    ylabel('\nabla c n')
                    title('Flux to bubble 1')
                    grid on
                    subplot(1,2,2)
                    plot(lRough(sum(PARAM.n(1:5))+1:end)-lRough(sum(PARAM.n(1:5))+1),conc(sum(PARAM.n(1:5))+1:end))
                    %legend('flux to bubble 2','Location','Best')
                    title('Flux to bubble 2')
                end
                xlabel('l')
                ylabel('\nabla c n')
                grid on

           end
           
           if plotConcLine==1
               
               %compute concentration
               Xline = linspace(-5,15,500);
               Yline = zeros(numel(Xline),1);
               [~,~,PHIfieldLine] = computeConcentrationField(Xline,Yline,x,y,conc,PARAM,0,coeffSub);
               
               if numel(manyTheta)>1 && numel(manyN)>1
                   if manyTheta(1)~=manyTheta(2) && manyXi(1)~=manyXi(2)
                        maxC(iii,kkk) = max(PHIfieldLine);
                   end
               elseif numel(manyN)>1
                   maxC(iii) = max(PHIfieldLine);
               elseif numel(manyTheta)>1
                   maxC(kkk) = max(PHIfieldLine);
               end
               
               %set concentration inside the bubble
               if numel(PARAM.n)==5
                   for k = 1:numel(PHIfieldLine)
                      if PHIfieldLine(k)==0
                          phiHenry = Hcc*(beta+2/rBubble);
                          PHIfieldLine(k) = phiHenry;
                      end
                   end
               elseif numel(PARAM.n)==6
                   for k = 1:numel(PHIfieldLine)
                      if PHIfieldLine(k)==0
                          if Xline(k)>xcmBubble-rBubble && Xline(k)<xcmBubble+rBubble
                            phiHenry = Hcc*(beta+2/rBubble);
                          else
                              phiHenry = Hcc*(beta+2/rBubble2);
                          end
                          PHIfieldLine(k) = phiHenry;
                      end
                   end
               end
               
               %plot
               if numel(manyN)==1 && numel(manyTheta)==1
                   figure
                   plot(Xline,PHIfieldLine)
                   hold on
                   grid on
                   xlabel('z')
                   ylabel('c')
                   title('Concentration along the axis')
               end
                
           end
           
           countPlot = countPlot+1;
    
end

end

%plot max conc
if isempty(maxC)==0
[~,ind] = min(abs(BC-maxC));
figure
if numel(manyTheta)>1 || numel(manyXi)>1
if manyTheta(1)~=manyTheta(2) && manyXi(1)==manyXi(2)
    plot(manyTheta,maxC)
    hold on
    plot(manyTheta(ind),maxC(ind),'.r','MarkerSize',40)
    xlabel('\theta')
    ylabel('max(c)')
    grid on
    title('max concentration along the axis')
elseif manyXi(1)~=manyXi(2) && manyTheta(1)==manyTheta(2)
    plot(manyXi,maxC)
    hold on
    plot(manyXi(ind),maxC(ind),'.r','MarkerSize',40)
    xlabel('\xi')
    ylabel('max(c)')
    grid on
    title('max concentration along the axis')
elseif manyXi(1)~=manyXi(2) && manyTheta(1)~=manyTheta(2)
    [~,h1] = contourf(manyTheta,manyXi,maxC,500);
    xlabel('\theta')
    ylabel('\xi')
    %grid on
    title('max concentration along the axis')
    colorbar
    set(h1,'LineColor','none')
end
end
end

%plot crirical radius
if plotCriticalRadius==1 && plotConcLine==1
    
    %crtical radius
    rCrit = 2*Hcc./(PHIfieldLine-beta*Hcc);
    
    figure
    semilogy(Xline,rCrit)
    xlabel('z')
    ylabel('r_{crit}')
    title('Critical radius of nucleation')
    grid on
    
    %flux to bubble
    r = 0.1:0.1:0.5;
    Cbubble = Hcc*(beta+2./r');
    
    %difference betweeen bubble and infinity
    %PHIfieldLine = conc(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3)))';
    %Xline = x{3};   Xline = (Xline(1:end-1)+Xline(2:end))/2;
    deltaC = repmat(Cbubble,1,numel(PHIfieldLine))-repmat(PHIfieldLine,numel(r),1);
    R0 = repmat(r',1,numel(PHIfieldLine));
    
    %flow to the bubble
    Qanalytical = -4*pi*deltaC.*R0;
    
    figure
    plot(Xline,Qanalytical)
    xlabel('z_0')
    ylabel('Q')
    grid on
    
    
end










