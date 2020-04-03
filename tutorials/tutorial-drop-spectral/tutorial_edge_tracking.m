%nonlinear BEM based solver for one droplet motion
%spectral with Chebyshev or Legendre

clear variables
close all

%% ADD LIBRARIES
REPOSITORY_NAME = '~/Documents/MATLAB/';    % path to the repository
addpath('../utils_one_drop_spectral')
add_paths_tutorials_drop_spectral(REPOSITORY_NAME);

%% PATH TO RESULTS
results = '../tutorial_results';
here = pwd;

tic

%% PHYSICAL PARAMETERS
PARAM.Ca = 0.1;                                       % capillary number
PARAM.visc = 5;     PARAM.Qdeflation = 0;           % viscosity ratio and deflation
V0 = 4/3*pi;                                        % drop volume
PARAM.placeShapeXCM = 0;                            % when doing edge tracking, replace always droplet in x0=0 at t0
PARAM.CaNL = 0;                                     % non linear extensional flow

%% BOUNDARY CONDITIONS
PARAM.BC = 1;                                       % BC: 1 is linear extensional flow, 2 is rising droplet, 3 is nonlinear extensional flow

%% INITIAL CONDITION
PARAM.uploadShape = 0;                                       % 0 the initial shape is analytical, 1 is from continuation, 2 the initial shape is from newton Method,  3 is from edge tracking, 4 is from DNS, 5 is from Minimal Seed
PARAM.D = [0.2, 0.6];                                              % deformation parameter of the initial condition
PARAM.shapeEllipse = 1;         % 0 initial shape is sphere plus P2 Legendre, 1 initial shape is an ellipse, 2 is pertubed with 2 legendre polynomials P2 P4, 3 is asymmetric shape by P3 Legendre, 4 is P2, P3

%% FRAME OF REFERENCE
PARAM.dropFrame = 1;                                % choose frame of reference: 0 is lab frame, 1 is drop frame
PARAM.Unormal = 1;                                  % 0 advect nodes with lagrangian velocity 1 with normal velocity, 3 use maesh stab, 4 radial velocity
PARAM.volume = 0;   PARAM.volTol = 1e-14;           % impose volume conservation with first mode

%% TIME DISCRETIZATION
PARAM.ODE = 2;          % 1 is ODE45, 2 is RK2, 3 is ODE23s, 4 is ODE23, 5 is ODE113, 6 is ODE23t, 7 is ODE15s, 8 is OD23tb
Tstart = 0;             % beginning of simulation
Tend = 500;             % end of simulation
maxDT = 1e-2;           % maximum time step if adaptive, otherwise simply time step
initialDT = maxDT;      % set initial time step
PARAM.Tend = Tend;      % end time

%% SPACE DISCRETIZATION
PARAM.legendre = 1;                                     % choose 0 Chebyshev, 1 Legendre, 2 is Legendre-Lobatto
PARAM.dealiasing = 50;  ratioGrid = 2;                  % dealiasing
PARAM.tresholdCoeffs = 1e-2;     PARAM.Rescale = 1;     % option for accuracy
PARAM.n = round(ratioGrid*PARAM.dealiasing)+1;          % number of grid points

%% REMESH
PARAM.remesh = 1;           % if 1 remesh is activated with newton method (spectral accuracy), 2 is with spline (lose spectral accuracy)
PARAM.remeshUniform = 0;    % impose uniform distribution of the mesh
PARAM.chooseMapping = 1;  PARAM.howStrong = 0.85; % 1 is normal mapping (like spectral parametrization), 2 is clustered toward the middle (how strong set how stronly to cluster)
PARAM.tolRemesh = 1e-5;     % tolerance for remeshing
PARAM.firstTimeRemesh = 0;  % don't perform remesh before this instant
PARAM.normRemesh = 1;       % norm on which I base the remeshing criteria: 1 is on the first derivative, 2 is on  the second derivative, 3 is on the sum of the two
PARAM.remeshSolve = 1;      % 1 newton method or 2 minimization fot remeshing
PARAM.MaxIter = 50;         % maximum iteration of "fsolve" or "fmincon"
PARAM.remeshStart = 0;      % remesh when drawning initial shape
PARAM.checkRemeshFail = 1;  % check if newton method for remeshing fails
PARAM.remeshMapping = chooseSpectralMapping(PARAM.chooseMapping,PARAM.howStrong);

%% SAVING OPTIONS
PARAM.saveData = 1;                                 % choose if to save
PARAM.saveDataIte = 0;                              % choose if to save every iteration
SaveHowMany = 100;                                  % save how many times
Tsave = linspace(Tstart,Tend,SaveHowMany+1);        % output at those time

%% EDGE TRACKING PARAMETERS
PARAM.edgeLoop = 20;                            % loop of edge tracking
PARAM.edgeStartingLoop = 1;                     % starting loop of edge tracking
PARAM.bisection = 2;                            % 1 average points resizing for volume conservation, 2 as before but resize also to bisect the area value
PARAM.deltaEdge = 1e-2;                         % next loop starts when a norm is smaller that this
PARAM.normEdge = 2;                             % on which norm I do the check: 1 is L-2 norm, 2 id L-infinity norm, 3 is excess area norm
PARAM.Dedge = [0.4 0.6];                        % deformation parameter first two shapes to check
PARAM.overlapModeEdge = [1 1];                  % kind of mode for perturbing first two shapes to check
PARAM.whichModeEdge = [2 2];                    % exactly which mode
PARAM.bruteForceBisection = 0;                  % ingnore bisection criterion, just bisect between two previous shapes
PARAM.deltaModify = 1;

%% CONVERGENCE
PARAM.checkRes = 1;              % check residuals or not (normal velocity to the interface)
PARAM.converge = 0;              % when resilduals are smaller, stop simulation
PARAM.convergeShape = 0.01;      % convergence based on shape (check convergence based on drolet shape, insert negative number if you want to deactivate)
PARAM.convergeShapeEdge = 0.02;  % to understand which trajectory is stable and which not

%% CHOOSE FILENAME
PARAM.algorithm = 2;    % edge tracking
PARAM.filename = chooseFilename(PARAM,maxDT,Tend,0);
PARAM.filenameIte = PARAM.filename;
PARAM.res = results;
PARAM.this_folder = here;
PARAM.V0 = V0;

%% INTEGRATION WEIGTHS AND DIFFERENTIATION MATRICES
if PARAM.legendre==1
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreIntDiff(0,1,PARAM.n+1);
    [PARAM.tcheb,PARAM.D1cheb,PARAM.D2cheb,PARAM.WGcheb] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==0
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.TTT] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==2
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreLobattoIntDiff(0,1,PARAM.n+1);
end

%% EDGE TRACKING
[Tedge,Yedge] = edgeTracking(Tsave,maxDT,initialDT,PARAM,V0);

simulationTime = toc;

%% SAVE RESULTS AND CLOSE
cd(results)
disp('Save data')
save(PARAM.filename)
cd(here)

%the end
disp('The End')















