%compute the motion of a conical motor with diffusion of chemical species

clear variables
close all

%% Add libraries paths
REPOSITORY_NAME = '~/Documents/MATLAB/';    % path to the repository
add_paths_tutorials(REPOSITORY_NAME);

%% Path to results
PARAM.res = '../tutorial_results';
PARAM.here = pwd;

%% Non dimensional parameters
L = 10;                         % channel lenght
H = 1;                          % chennel height
hFilm = 0.1;                    % Film thickness for the initial shape
PARAM.visc(4) = 10;             % viscosity ratio
PARAM.Bond = 0;                 % Bond number
PARAM.Ca = 0.05;                % Capillary number
alpha = 1.1;                    % Droplet size
PARAM.D(4) = 0;                 % Deformation parameter (make it ellpsoidal)
xStart0 = 0;                    % starting position of bubble
PARAM.Qsource = [0 0];          % impose bubble growth rate

%% Options
PARAM.cfunction = 0;                            % if 1, use c function for speed up (might be buggy)
PARAM.STstokes = 1;                             % if 1, use ad hoc method to integrate singular Green's functions
PARAM.kernelFreeSpace = 1;  PARAM.posWall = []; % if 1, use free space kernel, otherwise wall-bound kernel

%% Frame of reference
PARAM.dropFrame = 2;    % 0 is lab frame, 1 id drop frame, 2 in drop and replace center of mass in starting point

%% Time discretization
Tstart = 0;             % time at start
Tend = 500;             % time to end simulation
dt = 0.1;               % time step

%% Saving options
SaveHowMany = 1e2;                                  % save output how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time
PARAM.SaveDataIte = 1;                              % save dataat each checkpoint

%% Numerical discretization
nDrop = 40;                             % number of mesh elemnts for droplet
nPerLenght = 10;                        % number of mesh element on inlet/outlet/wall per unit lenght
PARAM.panels = [3 1];                   % number of panels per block (one block is a closed line, 1 panel is each boundary where BCs are applied)
PARAM.rotate = [0 0 0 0]/180*pi;        % geometry element rotation        
PARAM.n = [round(H*nPerLenght)
    round(L*nPerLenght)
    round(H*nPerLenght)
    nDrop];                             % number of element per panel
PARAM.orderVariableStokes = [0 0 0 1];    % 0 is constant on the elmennt, 1 is linear on the element
PARAM.orderGeometryStokes = [0 0 0 1];    % 0 is straight, 1 is curved (spline)
PARAM.SPlinesType = [0 0 0 2];            % 1 is natural splines, 2 is symmetric on the axis
PARAM.panelType = [0 0 0 2];              % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.panelInflate = [0 0 0 0];           % 0 is not inflating, 1 is inflating
PARAM.blockType = [0 2];                  % 0 is fixed wall, 1 is moving wall, 2 is droplet (this is for adding the force free equation)
PARAM.deflationBlock = [0 0];             % whether to apply deflation which is necessary for 0 viscoisty ratio and moving rigid objects (no by default)
PARAM.deflationConstant = [0 4*pi];       % constants associated to deflation

%% Geometry parameters
PARAM.geometryPanel = [0 0 0 2];        % 0 is a straight line, 1 ia an arc, 2 is a defromable object
PARAM.xStart = [-L/2 -L/2 L/2 nan];     % x starting point for the straight lines
PARAM.xEnd = [-L/2 L/2 L/2 nan];        % x ending point for the straight lines
PARAM.yStart = [0 H H nan];             % y starting point for the straight lines
PARAM.yEnd = [H H 0 nan];               % y ending point for the straight lines
PARAM.thetaStart = [nan nan nan 0];     % theta starting point for the arc
PARAM.thetaEnd = [nan nan nan pi];      % theta starting point for the arc
PARAM.rArc = [nan nan nan alpha];       % theta starting point for the arc
PARAM.x0_Circle = [nan nan nan xStart0];% center of the circle (axial direction)
PARAM.y0_Circle = [nan nan nan 0];      % center of the circle (radial direction)
PARAM.xCrotate = [0 0 0 0];             % rotation for circles
PARAM.yCrotate = [0 0 0 0];             % rotation for circles
PARAM.ellipseShape = [0 0 0 1];         % rotation for circles

% build geometry
[xInitial,yInitial,PARAM,tParametricBase] = buildGeometryPanelsGeneral(PARAM);
[xInitial{4},yInitial{4}] = drawBubbleInChannel(xStart0,0,alpha,hFilm,H,nDrop);

% vector of initial condition
initialXY = zeros(2*numel(xInitial{4}),1);
initialXY(1:2:end-1) = xInitial{4};
initialXY(2:2:end) = yInitial{4};

%% Boundary conditions
PARAM.typeBCstokes = [8 0 0 7];           % 1 is prescribed normal velocity
                                          % 2 is prescibed normal stress
                                          % 3 is prescribed tangent velocity
                                          % 4 is prescribed axial velocity
                                          % 5 is prescribed normal velocity and tangent stress
                                          % 6 is rigid body motion
                                          % 7 is prescribed normal stress due to curvature and gravity
                                          % 8 is prescribed axial stress and radial velocity
PARAM.addFlow = 0;                        % add background flow (for example extensional flow)
PARAM.velBCaxial = {nan 0 2*PARAM.Ca*(1-(yInitial{3}(1:end-1)+yInitial{3}(2:end))/2.^2) nan};
PARAM.velBCradial = {0 0 0 nan};
PARAM.stressBC = {0 nan nan 1};

%% Change deflation block when viscosity ratio is low
if PARAM.visc(4)<0.1
    PARAM.deflationBlock(2) = 1;
end

%% Repulsive forces
PARAM.repulsiveForces = [0 0];  % 0 is off, otherwise choose which repulsice force to use
PARAM.repulsiveOn = 2e-1;       % repulsive forces is on when distance is < than this
PARAM.coeffRepulsive = 1e-2;    % intensity of the repulsice forces
PARAM.smoothingRep = [];        % choose if to smooth out repulsice forces over the boundary

%% Remesh
PARAM.normRemesh = [1e-4 1e-4 1e-4 1e-4];           % When small, remesh more often
PARAM.remeshType = [0 1 0 4];                       % 1 element is split into two parts when coming close to another block, 2 uses a distribution, 4 remesh with a certain distributio (usually for deformable objects)
PARAM.remeshStep = 1;                               % minimum time integration step before repeating remesh
PARAM.remeshProximity = {4 4 4 2};                  % in case remesh by proximity, indicate which with panel to chek proximity
PARAM.maxElem = 1.1*[H/PARAM.n(1)
    L/PARAM.n(2)
    H/PARAM.n(3)
    alpha*pi/PARAM.n(4)];                           % max size of the elemnts (per panel)
PARAM.distActivateRemesh = [1 1 1 1];       
PARAM.adaptCoeff = [4 4 4 4];
PARAM.minSizeElemRemesh = [1e-3 1e-3 1e-3 1e-3]/2;  % min size of the elemnts (per panel)
PARAM.coeffDist = 1;                                % how much smaller than distnce it has to be
PARAM.distr = [nan nan nan 2];                      % choose distr, 1 is uniform, 2 cluster more close to upper wall
PARAM.adaptDistr = [nan nan nan 1];                 % how strong to use distribution
PARAM.maxNumberTotalElem = 1e3;                     % total max number of elements admitted

%% Print info to screen
printToScreenStokes(PARAM)              % info related to Stokes equation (e.g BCs)
printToScreenTimeStepping(0,dt,0)       % info time stepping
printToScreenRemesh(PARAM)              % info remesh
printToScreenVerticalTube(L,PARAM.x0_Circle(4),alpha,PARAM.Bond,PARAM.visc(4),PARAM.dropFrame,PARAM.Ca) % infor on current non-dimensionalization

%% Create filename for saving results
PARAM.filename = chooseFilenameVerticalTube(PARAM,sum(PARAM.n),L,Tend,dt,alpha,PARAM.x0_Circle(4),PARAM.Bond,0,PARAM.visc(4),sum(PARAM.repulsiveForces),PARAM.Ca);

%% Handle functions
fRising = @(t,var) computeVelocityDropInPipe(t,var,tParametricBase,PARAM);  % MOST IMPORTANT: Compute the interface velocity as a function it it position dx/dt=F(x)
remeshFunction = @(t,var) remeshPanelsOneBubble(t,var,PARAM);           % Remesh handle function
volCorrFunction = @(t,var) volCorrBubble(t,var,4/3*pi*alpha^3);         % Volume correction handle function
    
%% Saving and event function
outFun = @(t,var,T,Y,V) eventBlocksRisingInTubeRK2(t,var,T,Y,V,tParametricBase,PARAM);

tic

%% Run time stepping RK2
[T,Y,V] = RK2mostGeneralBubbleVolCorr(fRising,Tsave,initialXY,dt,dt,outFun,remeshFunction,volCorrFunction,PARAM.dropFrame==2);

simulationTime = toc;

%% Clear function handles, save results and close
clear outfun
clear fRising
clear remeshFunction
clear volCorrFunction

%save results
disp('Save results')
cd(PARAM.res)
save(PARAM.filename)
cd(PARAM.here)

disp('The end')


