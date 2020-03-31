%post processing of a cone moving because of phoretic effects

clear variables
close all

%parameters
manyN = 200*ones(1,4);
manyNmodes = [3 4 5 10];
M = 1;
res = 1;

%options
range = 3;  shift = 0;
plotFieldConc = 1;
plotSutfaceConc = 1;
plotVelField = 0;   labFrame = 1;   % 0 is in motor frame, 1 is in lab frame
checkDecay = 1;     xDecay = logspace(1,5,100); yDecay = zeros(1,numel(xDecay));

%grid for field quantities
nGrid = 20;
x = linspace(-range+shift,range+shift,2*range*nGrid);
y = linspace(0,range,range*nGrid);
[X,Y] = meshgrid(x,y);

%upload directory
if res==1
   dir = '~/Documents/MATLAB/phoreticSwimmer/results/'; 
end

%initialize
U0 = zeros(numel(manyN),1);
U0plot = zeros(numel(manyN),1);
cellLegend = cell(numel(manyN),1);

countPlot = 1;
doNotPlot = 0;
for i = 1:numel(manyN)
    
    display([num2str(i) ' of ' num2str(numel(manyN))])
    
   %parameters
   n = manyN(i);
   nModes = manyNmodes(i);
   
   %upload data
   pathUp = [dir 'spherePhoretic_n=' num2str(n) '_nMode=' num2str(nModes) '.mat'];
   load(pathUp)
   
   U0(i) = U;
   
   %plot motor
   if checkDecay==1

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
                cellLegend{countPlot} = ['nModes=' num2str(manyNmodes(i))];

   end
            
   figure
   if plotFieldConc==1
                
               %compute concentration
               [X,Y,PHIfield] = computeConcentrationField(X,Y,x,y,conc,PARAM);
               contourf(X,Y,PHIfield)
               hold on
               contourf(X,-Y,PHIfield)
               colorbar
                
   end

   frame = '';
   if plotVelField==1
                
   %compute concentration
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
                
   end
           
   %plot cone
   plotGeometry(x,y,PARAM);
   grid on
   axis equal
   axis([-range+shift range+shift -range range])
   title([frame ' U=' num2str(yStokes(end)) ' ' num2str(nMode) ' modes'])

   if plotSutfaceConc==1
       
        %compute location of the singularity
        PARAM.orderVariable = PARAM.orderVariableLaplace;
        PARAM.orderGeometry = PARAM.orderGeometryLaplace;
        [Xsing,Ysing] = computeSingularityLocation(x,y,PARAM);

        figure
        lRough = [0; cumsum(sqrt(diff(Xsing).^2+diff(Ysing).^2))];
        plot(lRough,conc)
        hold on
        plot(lRough,vSlip)
        xlabel('l')
        ylabel('c,\nabla c')
        grid on
        title(['U=' num2str(yStokes(end)) ' ' num2str(nMode) ' modes'])

   end
           
   U0plot(countPlot) = U0(i);
   countPlot = countPlot+1;
    
end

%plot swimming velocity
figure
plot(manyNmodes,abs(U0),'-o')
hold on
xlabel('modes')
ylabel('|U|')
title('Motor Velocity')
grid on

if checkDecay==1
    figure(1)
    legend(cellLegend,'Location','Best')
end











