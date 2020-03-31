%post processing of test time stepping spectral

close all
clear variables

%data
CaUP = [];
BondUP = 5;
viscUP = 1;
n = 50;
TendUP = 20;
ODEup = 0;          % 1 id ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
BC = 7;             % 1 is extensional flow, 2 is rising droplet
dtUP = 1e-3;
IDup = 1:2000;         % initial shape
res = 0;            % results, server or results/forThesis
numCoeffIC = 2;
nameTag = ['RandomSearch_' num2str(numCoeffIC) 'modes_'];
extractFromTransient = 2;

% parametric study on this
param = IDup;
%param = 38;

%option
plotShape = 0;
plotElongNorm = 0;
zoom = 0;
plotRes = 0;
plotAreaNorm = 0;
Tbreak = 2*TendUP;
plotVolErr = 0;
plotPhaseSpace = 1;

%directory
if res==0
    dir = '~/Documents/MATLAB/droplet_simulations/server/';
elseif res==1
    dir = '~/Documents/MATLAB/droplet_simulations/results/';
elseif res==2
    dir = '~/Documents/MATLAB/droplet_simulations/results/forThesis/manyExtensionalFlow/';
end

%initialize
breakOrNot = cell(2,1);

for k = 1:numel(param)
    
    display([num2str(k) ' of ' num2str(numel(param)) ', ID=' num2str(param(k))])
    
    Dup = param(k);

    %filename
    if isempty(BondUP)
        name = [nameTag 'oneDropBEM_ODE=' num2str(ODEup) '_n=' num2str(n) '_BC=' num2str(BC) '_Ca=' num2str(CaUP) '_visc=' num2str(viscUP) '_D=' num2str(Dup) '_maxDT=' num2str(dtUP) '_Tend=' num2str(TendUP) '.mat'];
    elseif isempty(CaUP)
        name = [nameTag 'oneDropBEM_ODE=' num2str(ODEup) '_n=' num2str(n) '_BC=' num2str(BC) '_Bo=' num2str(BondUP) '_visc=' num2str(viscUP) '_ID=' num2str(Dup) '_maxDT=' num2str(dtUP) '_Tend=' num2str(TendUP) '.mat'];
    end

    %upload data
    upload = [dir name];
    try
    load(upload)

    T = allRes{1};
    YYY = allRes{2};
    VVV = allRes{3};

    %number of iteration
    ite = round(numel(T));

    %initialize
    normElon = zeros(ite,1);
    elemNum = zeros(ite,1);
    res = zeros(ite,1);
    V = zeros(ite,1);
    Area = zeros(ite,1);
    Evisc = zeros(ite,1);
    Dellipse = zeros(ite,1);
    xcm = zeros(ite,1);
    Anorm = zeros(ite,1);
    f1 = zeros(ite,1);
    f2 = zeros(ite,1);
    f3 = zeros(ite,1);
    V0 = 4/3*pi;

    %loops
    plotCount = 1;
    shift = 0;  islast = 0;
    fFromTransient = zeros(ite-1,extractFromTransient);
    for i = 1:ite

       %get currents modes
       Y = YYY{i};
       x{1} = Y(1:2:end-1)';
       y{1} = Y(2:2:end)';
       xGrid = x{1};
       yGrid = y{1};
       
       if plotShape==1
            figure(10) 
            plotGeometryStokes(x,y,0,[],[],[],0,PARAM)
            xyLabelTex('z','r')
            hold off
            title('Drop shape')
      
            axis equal
            axis([-4 2 -2 2])
            if BC==7
                camroll(90)
            end
            drawnow
       end


       elemNum(i) = numel(xGrid);

       %compute elongation norm
       xcm(i) = centerOfMassBlockAxis(x,y,1,PARAM);
       RforNorm = sqrt((xGrid-xcm(i)).^2+yGrid.^2);
       normElon(i) = abs(RforNorm(end)-xcm(i));

       %compute volume abd surface area
       V(i) = axis_int_gauss_vect(xGrid,yGrid);
       Area(i) = surf_gauss_vect(xGrid,yGrid);
       
       if i==1
           xHere = xGrid-xcm(i);
           %xHere = xGrid;
           radius = sqrt(xHere.^2+yGrid.^2);
           theta = atan(yGrid./xHere);
           theta = theta + pi*(theta<0);
           fHere = LegendreSeriePolar(theta',radius',max(numCoeffIC,extractFromTransient)+1,0);
           if min(yGrid)<0
               clear fFirst
               break
           end
           fFirst = fHere(2:end);
           %fFirst = PARAM.f{1};
       end
       
       %compute modes of transient (polar coordinates)
       if i>1 && extractFromTransient>0
           xHere = xGrid-xcm(i);
           %xHere = xGrid;
           radius = sqrt(xHere.^2+yGrid.^2);
           theta = atan(yGrid./xHere);
           theta = theta + pi*(theta<0);
           fHere = LegendreSeriePolar(theta',radius',extractFromTransient+1,0);
           fFromTransient(i-1,:) = fHere(2:end);
       end

       if plotAreaNorm==1

           R0 = nthroot(V(i)/4/pi*3,3);
           Anorm(i) = (Area(i) - 4*pi*R0^2)./(Area(1) - 4*pi*R0^2);

       end

       %defromation parameter
       L = max(xGrid)-min(xGrid);
       B = 2*yGrid(round((PARAM.n+1)/2));
       Dellipse(i) = (L-B)/(L+B);

       %get time
       t = T;

       if PARAM.ODE==0

           Vhere = VVV{i};
           UxN = Vhere(1:2:end-1);
           UyN = Vhere(2:2:end);
           UnABS = sqrt(UxN.^2+UyN.^2);
           res(i) = norm(UnABS,Inf);

       else
           error('Not implemented')
       end

       %stop when there are no more data
       if i~=ite
       if T(i+1)>Tbreak||(T(i+1)==0&&k>1)
           V = V(1:i);
           t = T(1:i);
           res = res(1:i);
           normElon = normElon(1:i);
           disp('Break')
           
           breakOrNot{1} = [breakOrNot{1}; fFirst];
           if extractFromTransient>0
               [rowF, ~] = size(fFromTransient);
               breakOrNot{1} = [breakOrNot{1}; fFromTransient];
           end
           if abs(Dellipse(i))<1e-2
                breakOrNot{2} = [breakOrNot{2}; 0];
                if extractFromTransient>0
                    breakOrNot{2} = [breakOrNot{2}; zeros(rowF,1)];
                end
           else
                breakOrNot{2} = [breakOrNot{2}; 1];
                if extractFromTransient>0
                    breakOrNot{2} = [breakOrNot{2}; ones(rowF,1)];
                end
           end
           
           break;
       end
       end

    end

    if plotVolErr==1

        errV = (V-V0)/V0;
        figure(1)
        plot(t,errV)
        hold on
        xlabel('t')
        ylabel('err_V')
        grid on
        title('Error On Volume')

    end

    if plotRes==1
        
        figure(2)
        semilogy(t,res)
        hold on
        %xlabel('t')
        %ylabel('res')
        xyLabelTex('t','||\mathbf{u} \cdot \mathbf{n}||_\infty')
        title('residuals')
        grid on

    end

    if plotElongNorm==1

        figure(3)
        %semilogy(t,normElon,'k')
        plot(t,normElon-1)
        hold on
        %xlabel('t')
        %ylabel('L/a')
        xyLabelTex('t','L')
        title('Stone elongation')
        grid on
        %axis([0 200 1 10])

    end

    if plotAreaNorm==1

        figure(4)
        semilogy(t(1:i),Anorm(1:i))
        hold on
        %xlabel('t')
        %ylabel('\Delta S')
        xyLabelTex('t','\Delta S')
        title('Excess area')
        grid on

    end
    
    breakOrNot{1} = [breakOrNot{1}; fFirst];
    if extractFromTransient>0
        [rowF, ~] = size(fFromTransient);
        breakOrNot{1} = [breakOrNot{1}; fFromTransient];
    end
    if abs(Dellipse(i))<1e-2
        breakOrNot{2} = [breakOrNot{2}; 0];
        if extractFromTransient>0
             breakOrNot{2} = [breakOrNot{2}; zeros(rowF,1)];
        end
    else
        breakOrNot{2} = [breakOrNot{2}; 1];
        if extractFromTransient>0
             breakOrNot{2} = [breakOrNot{2}; ones(rowF,1)];
        end
    end

    %print simulation time
    %disp(['Simulation time T=' num2str(simulationTime/60) ' minutes'])
    
    catch
        disp('No data');
    end
    
end

%plot phase space, break or not
if plotPhaseSpace==1
    
    sizeF = size(breakOrNot{1});
    if sizeF(2)==2
        
        fff = breakOrNot{1};
        tag = breakOrNot{2};
        f1break = fff(tag==1,1);
        f1noBreak = fff(tag==0,1);
        f2break = fff(tag==1,2);
        f2noBreak = fff(tag==0,2);
        
        figure
        plot(f1break,f2break,'x')
        hold on
        plot(f1noBreak,f2noBreak,'o')
        xyLabelTex('f_1','f_2')
        grid on
        legend('Breakup','No breakup','Location','Best')
        %axis([0 1 0 1])
        
    elseif sizeF(2)==3
        
        fff = breakOrNot{1};
        tag = breakOrNot{2};
        f1break = fff(tag==1,1);
        f1noBreak = fff(tag==0,1);
        f2break = fff(tag==1,2);
        f2noBreak = fff(tag==0,2);
        f3break = fff(tag==1,3);
        f3noBreak = fff(tag==0,3);
        
        figure
        plot3(f1break,f2break,f3break,'x')
        hold on
        plot3(f1noBreak,f2noBreak,f3noBreak,'o')
        xyLabelTex3d('f_1','f_2','f_2')
        grid on
        legend('Breakup','No breakup','Location','Best')
        %axis([0 1 0 1])
        
    end

end









