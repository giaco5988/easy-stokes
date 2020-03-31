%post processing of a cone moving because of phoretic effects

clear variables
close all

%parameters
nTheta = 100;
manyN = 100*ones(1,nTheta);
manyXi = 1;
manyHpost = 0.05*ones(1,nTheta);
%manyRpost = 0.05:0.01:1;
manyRpost = 0.25;
manyThetaPost = linspace(1e-3,pi-1e-3,nTheta);
%manyTheta = 2.5;
flux = [0 0 -1 0];
M = 0;
res = 1;
parametric = 0;
rec = 0;

%options
range = 2;  shift = -0.5;
plotSnaphot = 1;
if numel(manyThetaPost)>1
    %thetaPlot = [manyTheta(1) manyTheta(6) manyTheta(26) manyTheta(51) manyTheta(79) manyTheta(end)];
    %thetaPlot = [manyTheta(1) manyTheta(17) manyTheta(50) manyTheta(80)];
    thetaPlot = [manyThetaPost(2) manyThetaPost(5) manyThetaPost(10) manyThetaPost(15) manyThetaPost(20)];
    %thetaPlot = manyThetaPost([5:5:100]);
elseif numel(manyThetaPost)==1
    thetaPlot = manyThetaPost(1);
end
plotConc = 0;
plotSutfaceConc = 0;
plotVelField = 0;   labFrame = 1;   % 0 is in motor frame, 1 is in lab frame
checkDecay = 0;     xDecay = logspace(1,5,100); yDecay = zeros(1,numel(xDecay));
plotVelAlongAxis = 1;

%grid for field quantities
nGrid = 10;
x = linspace(-range+shift,range+shift,2*range*nGrid);
y = linspace(0,range,range*nGrid);
[X,Y] = meshgrid(x,y);

%upload directory
if res==1
   dir = '~/Documents/MATLAB/phoreticSwimmer/results/'; 
elseif res==0
   dir = '~/Documents/MATLAB/phoreticSwimmer/server/';
end

if parametric==1
    if numel(manyRpost)==1 && manyRpost~=0
        pathUp = [dir 'conePhoreticParametric_n=' num2str(manyN(1)) '_L=' num2str(manyXi(1)) '_r=' num2str(manyRpost) '_h=' num2str(manyHpost(1)) '_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_f4=' num2str(flux(4)) '.mat'];
    elseif numel(manyRpost)==1 && manyRpost==0
        pathUp = [dir 'closedConeParametric_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_n=' num2str(manyN(1)) '_xi=' num2str(manyXi(1)) '_h=' num2str(manyHpost(1)) '.mat'];
    elseif numel(manyRpost)>1
        pathUp = [dir 'conePhoreticParametric_n=' num2str(manyN(1)) '_L=' num2str(manyXi(1)) '_h=' num2str(manyHpost(1)) '_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_f4=' num2str(flux(4)) '.mat'];
    end
    load(pathUp)
end

%initialize
manyU0 = zeros(numel(manyThetaPost),numel(manyRpost));
manySnum = zeros(numel(manyThetaPost),numel(manyRpost));
U0plot = zeros(numel(thetaPlot),1);
Splot = zeros(numel(thetaPlot),1);
cellLegend = cell(numel(thetaPlot),1);

if manyRpost==0
    closed = 1;
else
    closed = 0;
end

for kkk = 1:numel(manyRpost)
    
countPlot = 1;
doNotPlot = 0;
r = manyRpost(kkk);
n = manyN(kkk);

if numel(manyRpost)>1
    display([num2str(kkk) ' of ' num2str(numel(manyRpost))])
end

for iii = 1:numel(manyThetaPost)
    
    if numel(manyRpost)==1
        display([num2str(iii) ' of ' num2str(numel(manyThetaPost))])
    end
    
   %parameters
   h = manyHpost(iii);
   
   %upload data
   if M==0 && closed==0 && parametric==0 && rec==0
        pathUp = [dir 'conePhoreticParametric_n=' num2str(manyN(1)) '_L=' num2str(manyXi(1)) '_r=' num2str(manyRpost(1)) '_h=' num2str(h) '_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_f4=' num2str(flux(4)) '.mat'];
        load(pathUp)
   elseif M==1 && closed==0 && parametric==0 && rec==0
        pathUp = [dir 'conePhoreticAllMob_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_f4=' num2str(flux(4)) '_n=' num2str(n) '_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
        load(pathUp)
   elseif M==0 && closed==1 && parametric==0 && rec==0
        pathUp = [dir 'closedCone_f1=' num2str(flux(1)) '_f2=' num2str(flux(2)) '_f3=' num2str(flux(3)) '_n=' num2str(n) '_xi=' num2str(manyXi) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
        load(pathUp)
   elseif rec==1 && parametric==0
        pathUp = [dir 'coneReciprocalTheorem_n=' num2str(n) '_xi=' num2str(manyXi(1)) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];
        load(pathUp)
   end
   
   theta = manyThetaPost(iii);
   
   %solution
   yLaplace = solutionLaplace{iii,kkk};
   yStokes = solutionStokes{iii,kkk};
   x = xxx{iii,kkk};
   y = yyy{iii,kkk};
   
   if parametric==0 && rec==0
       %motor velocity
       manyU0(iii,kkk) = yStokes(end);
       %compute stresslet
       manySnum(iii,kkk) = computeStressLetAxisNumerical(x,y,yStokes(1:end-1),PARAM);
   elseif parametric==1 && rec==0
       %motor velocity
       manyU0(iii,kkk) = U0(iii,kkk);
       %compute stresslet
       manySnum(iii,kkk) = Snum(iii,kkk);
   end
   %plot motor
   if plotSnaphot==1 && doNotPlot==0
       
       for qqq = 1:numel(PARAM.n)
                    PARAM.n(qqq) = numel(x{qqq})-1;
       end
       
       if thetaPlot(countPlot)==manyTheta(iii)
           
           if plotVelAlongAxis==1
       
              %compute velocity along the axis
              Xaxis = linspace(-4,4,400);
              Yaxis = zeros(1,numel(Xaxis));
              [~,~,uField,vField] = computeVelPressField(Xaxis,Yaxis,x,y,yStokes(1:end-1),zeros(2*numel(Xaxis),1),PARAM,yStokes(end),0);

              figure(1)
              hold on
              %Uabs = sqrt(uField.^2+vField.^2);
              plot(Xaxis,uField)
              xlabel('x')
              ylabel('U')
              grid on
              title('Velocity along the axis')
              drawnow
              
           end

            if checkDecay==1

                figure(1)
                if countPlot>1
                    hold on
                end
                if rec==0
                    [~,~,uField,vField] = computeVelPressField(xDecay,yDecay,x,y,yStokes(1:end-1),0,PARAM,yStokes(end),0);
                elseif rec==1
                    [~,~,uField,vField] = computeVelPressField(xDecay,yDecay,x,y,yStokes,0,PARAM,yStokes(end),0);
                end
                Uabs = sqrt(uField.^2+vField.^2);
                loglog(xDecay,Uabs,'o-')
                xlabel('d')
                ylabel('|U|')
                title('Front velocity decay')
                grid on
                cellLegend{countPlot} = ['\theta=' num2str(manyTheta(iii))];

            end
            
           figure
           if plotConc==1
                
               %compute concentration
               [X,Y,PHIfield] = computeConcentrationField(X,Y,x,y,yLaplace,PARAM);
               contourf(X,Y,PHIfield)
               hold on
               contourf(X,-Y,PHIfield)
               colorbar
                
           end

           frame = '';
           if plotVelField==1
                
               %compute velocity field
               if rec==0
                    [X,Y,uField,vField] = computeVelPressField(X,Y,x,y,yStokes(1:end-1),zeros(2*numel(X),1),PARAM,yStokes(end),1);
               elseif rec==1
                   [X,Y,uField,vField] = computeVelPressField(X,Y,x,y,yStokes,[PARAM.Uunder*ones(numel(X),1); zeros(numel(x),1)],PARAM,0,1);
               end
               
               if labFrame==0
                    uField = uField-yStokes(end);
                    frame = 'Motor frame';
               elseif labFrame==1
                    frame = 'Lab frame';
               end
               streamslice(X,Y,uField,vField)
               hold on
               uAbs = sqrt(uField.^2+vField.^2);
               contourf(X,-Y,uAbs,linspace(0,2,10))
               colorbar
               drawnow
                
           end
           
           %plot cone
           plotGeometry(x,y,PARAM);
           grid on
           axis equal
           xlabel('x')
           ylabel('r')
           axis([-range+shift range+shift -range range])
           if rec==0
                title([frame ' \theta=' num2str(theta) ' U=' num2str(yStokes(end))])
           elseif rec==1
                title(['\theta=' num2str(theta)])
           end

           if plotSutfaceConc==1

                figure
                lRough = [0; cumsum(sqrt(diff(Xsing).^2+diff(Ysing).^2))];
                plot(lRough,yLaplace)
                hold on
                for qqq = 1:numel(PARAM.n)
                    %finite differences
             
                    D1 = finiteDifference1D(PARAM.n(qqq),[2 0],1);
                
                    %get range
                    if qqq==1
                        rangeHere = 1:PARAM.n(1);
                    else
                        rangeHere = sum(PARAM.n(1:qqq-1))+1:sum(PARAM.n(1:qqq));
                    end
                    
                    if sum(PARAM.velBC{qqq})~=0
                        col = 'r';
                    else
                        col = 'k';
                    end
                    %compute slip velocity
                    vSlip = computeSlipVelPhoretic(x{qqq},y{qqq},yLaplace(rangeHere),PARAM,qqq,D1);
                    plot(lRough(rangeHere),vSlip,col)
                end
                xlabel('l')
                ylabel('c,\nabla c')
                grid on
                title(['\theta=' num2str(theta) ' U=' num2str(yStokes(end))])

           end
           
           U0plot(countPlot) = manyU0(iii,kkk);
           Splot(countPlot) = manySnum(iii,kkk);
           if countPlot==numel(thetaPlot)
                doNotPlot = 1;
           end
           countPlot = countPlot+1;
                
       end

   end
    
end

end

if numel(manyRpost)==1 && numel(manyTheta)>1

%plot swimming velocity
figure
if manyTheta(1)~=manyTheta(2)
    par = manyTheta;
    parAxis = '\theta';
elseif manyH(1)~=manyH(2)
    par = manyH;
elseif manyXi(1)~=manyXi(2)
    par = manyXi;
end
plot(par,manyU0,'-')
hold on
xlabel(parAxis)
ylabel('U')
title('Motor Velocity')
grid on
if plotSnaphot==1
   plot(thetaPlot,U0plot,'.','MarkerSize',30) 
end

if checkDecay==1
    figure(1)
    legend(cellLegend,'Location','Best')
end

%plot stresslet
figure
if manyTheta(1)~=manyTheta(2)
    par = manyTheta;
    parAxis = '\theta';
elseif manyH(1)~=manyH(2)
    par = manyH;
elseif manyXi(1)~=manyXi(2)
    par = manyXi;
end
plot(par,manySnum,'-')
hold on
xlabel(parAxis)
ylabel('S/A')
title('Stresslet intensity')
grid on
if plotSnaphot==1
   plot(thetaPlot,Splot,'.','MarkerSize',30) 
end

if checkDecay==1
    figure(1)
    legend(cellLegend,'Location','Best')
end

elseif numel(manyRpost)>1 && numel(manyTheta)>1
    
    %[X,Y] = meshgrid(manyTheta,manyRpost/L);
    
    figure
    contourf(manyTheta,manyRpost/L,manyU0')
    %contourf(X,Y,manyU0')
    xlabel('\theta')
    ylabel('h/L')
    title('Motor Velocity')
    colorbar
    
    figure
    contourf(manyTheta,manyRpost/L,manySnum')
    xlabel('\theta')
    ylabel('h/L')
    title('Stresslet intensity')
    colorbar
end











