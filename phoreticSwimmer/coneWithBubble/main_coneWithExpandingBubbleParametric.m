%compute the diffusion of chemical species around a Janus particle with BEM

clear variables
close all

%results
PARAM.res = '~/Documents/MATLAB/phoreticSwimmer/results';
PARAM.here = pwd;

%physical parameters
r = 1;                             % average cone radius
%manyRbubble = 0.01:0.01:0.05;        % radius of the bubble
manyRbubble = 0.01:0.005:0.1;
L = 10;                            % motor lenght
h = 0.2;                           % wall thickness
theta = 5/180*pi;                         % inclinitation of the cone
flux = [0 0 -1 0 0];               % BC for flux
manyInPos = 0:2.5:10;              % initial bubble position
Hcc = 0.1;
beta = 1;
Rcut = r;

%time marching
loop = 1000;
dt = 1e-3;

%options
PARAM.cfunction = 0;
PARAM.STlaplace = 1;
PARAM.STstokes = 1;

%geometry parameters
nPerLenght = 20;
PARAM.panels = [4 1];                              % panels per block
PARAM.geometryPanel = [0 1 0 1 1];                % 0 is a straight line, 1 ia an arc
PARAM.xStart = [L nan 0 nan nan];             % x starting point for the straight lines
PARAM.xEnd = [0 nan L nan nan];               % x ending point for the straight lines
PARAM.yStart = [r+h/2 nan r-h/2 nan nan];             % y starting point for the straight lines
PARAM.yEnd = [r+h/2 nan r-h/2 nan nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan pi/2 nan -pi/2 0];               % theta starting point for the arc
PARAM.thetaEnd = [nan 3*pi/2 nan pi/2 pi];               % theta starting point for the arc               % theta starting point for the arc
PARAM.y0_Circle = [nan r nan r 0];

%function profile for BCs
PARAM.fluxBC{1} = 0;
PARAM.fluxBC{2} = 0;
PARAM.fluxBC{3} = flux(3);
PARAM.fluxBC{4} = 0;
PARAM.fluxBC{5} = 0;
PARAM.concBC{1} = 0;
PARAM.concBC{2} = 0;
PARAM.concBC{3} = 0;
PARAM.concBC{4} = 0;

%numerics parameters for Stokes
PARAM.typeBCstokes = [3 3 3 3];           % 1 is prescribed normal velocity, 2 is prescibed normal stress, 3 is prescribed tangent velocity
PARAM.orderVariableStokes = [0 0 0 0];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 0];    % 0 is straight, 1 is curved (spline)
PARAM.panelType = [0 0 0 0];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.blockType = 1;                      % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)

%numerics Laplace
PARAM.typeBClaplace = [2 2 2 2 1];            % 1 is prescribed concentration, 2 is prescibed flux
PARAM.orderVariableLaplace = [0 0 0 0 0];     % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryLaplace = [0 0 0 0 0];     % 0 is straight, 1 is curved (spline)

%initialize
Tscale = zeros(numel(manyInPos),numel(manyRbubble));
manyQstart = zeros(numel(manyInPos),numel(manyRbubble));

%filename
PARAM.filename = ['coneExpandBubbleParametric_xi=' num2str(L) '_theta=' num2str(theta) '_h=' num2str(h) '.mat'];

for l = 1:numel(manyInPos)
    
    InPos = manyInPos(l);
    
    for k = 1:numel(manyRbubble)
        
    rBubble = manyRbubble(k);
    PARAM.x0_Circle = [nan 0 nan L InPos];
    PARAM.rArc = [nan h/2 nan h/2 rBubble];
    PARAM.n = [round(L*nPerLenght) round(h*pi/2*nPerLenght) round(L*nPerLenght) round(h*pi/2*nPerLenght) round(pi*rBubble*nPerLenght)+20];                  % number of element per panel
    
    % first iteration
    PARAM.concBC{5} = (2/rBubble + beta)*Hcc;    

    %build geometry
    xC = 0;   yC = r;
    tParametric = parametricGrid(PARAM);

    %build physical shape
    [x,y] = buildGeometryPanelsParametric(tParametric,PARAM);

    %rigid body rotation
    [x,y,tParametric,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd] = rigidBodyRotation(x,y,tParametric,xC,yC,theta,PARAM.x0_Circle,PARAM.y0_Circle,PARAM.xStart,PARAM.xEnd,PARAM.yStart,PARAM.yEnd,PARAM.thetaStart,PARAM.thetaEnd,PARAM.geometryPanel,1:4);

    %print to screen
    %printToScreenLaplace(PARAM)
    %printToScreenStokes(PARAM)

    %finite differences
    D1{1} = finiteDifference1D(PARAM.n(3),[2 0],1);

        for i = 1:loop

            display(['T=' num2str(dt*(i-1))])

            %solve laplace equation
            [yLaplace,Xsing,Ysing] = BEM_Laplace(x,y,PARAM);

            %compute slip
            [vSlip,larc] = computeSlipVelPhoretic(x{3},y{3},yLaplace(sum(PARAM.n(1:2))+1:sum(PARAM.n(1:3))),PARAM,3,D1{1});

        %     %plot flux
        %     figure
        %     plot(Xsing(sum(PARAM.n(1:4))+1:sum(PARAM.n)),yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n)))
        %     xlabel('x')
        %     ylabel('\nabla c')
        %     grid on
        %     title('Flux trough the bubble')

            %compute flow rate
            Q = computeFlowRate(x{5},y{5},yLaplace(sum(PARAM.n(1:4))+1:sum(PARAM.n)),5,PARAM);
            manyQ(i) = Q;
            manyV(i) =  axis_int_gauss_vect(x{5},y{5});
            
            if i==1
                
               manyQin(l,k) = Q;
                
            end

%             figure(10)
%             %plotGeometry(x,y,PARAM)
%             plot(x{5},y{5})
%             axis([-1 11 -2 2])
%             axis equal
%             hold off
%             drawnow
            
            %inflate bubble
            Un = Q/(3*beta^rBubble^2+4*rBubble);
            rBubble = rBubble + dt*Un;

            %new concentration on the bubble
            PARAM.concBC{5} = Hcc*(beta + 2/rBubble);

            %function profile for BCs
            PARAM.velBC{1} = 0;
            PARAM.velBC{2} = 0;
            PARAM.velBC{3} = vSlip;
            PARAM.velBC{4} = 0;

            %solve Stokes equation
            %[yStokes,Xsing,Ysing] = BEM_Stokes(x,y,PARAM);

            %new bubble geometry
            %PARAM.n(5) = round(10*rBubble*pi/2*nPerLenght);
            [x{5},y{5}] = buildArc(rBubble,InPos,0,[0 pi],PARAM,5);

            %bubble disapperas
            if min(y{5}<0.02)
                display('Bubble has vanished')
                vanishes = 1;
                Tscale(l,k) = i*dt;
                break;
            end

            %bubble is larger than the tube
            phi = FigInOut(x{5},y{5},[x{1} x{2} x{3} x{4}]',[y{1} y{2} y{3} y{4}]');
            comp = max(phi>pi);
            if comp || max(y{5})>Rcut
                display(['Bubble touches the tube of its radius is largger that r=' num2str(Rcut)])
                vanishes = 0;
                Tscale(l,k) = i*dt;
                break;
            end

        end
        
        if loop>1
            if vanishes==1
                    col = 'xk';
            elseif vanishes==0
                    col = 'or';
            end
        end
    
        if loop>1
        figure(3)
        hold on
        plot(manyInPos(l),manyRbubble(k),col)
        xlabel('Initial Position')
        ylabel('Initial bubble radius')
        title(['Bubble grows or disappears \theta=' num2str(theta)])
        grid on
        drawnow
        end

    end
    
end

%plot geometry
% figure
% hold on
% %plotGeometry(x,y,PARAM)
% axis equal
% axis(0.1*[-2 12 -3 3])
% grid on
% title('Motor Geometry')
% xlabel('z')
% ylabel('r')

if loop>1
    %plot
    figure
    hold on
    plot(manyInPos,Tscale)
    grid on
    title('Tscale')
    xlabel('z_0')
    ylabel('T')
end

%plot
figure
hold on
plot(manyInPos,manyQin)
grid on
title('Inflation rate')
xlabel('z_0')
ylabel('Q')

cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

display('The end')










