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
PARAM.D = 0.2;                                              % deformation parameter of the initial condition

%% FRAME OF REFERENCE
PARAM.dropFrame = 1;                                % choose frame of reference: 0 is lab frame, 1 is drop frame
PARAM.Unormal = 1;                                  % 0 advect nodes with lagrangian velocity 1 with normal velocity, 3 use maesh stab, 4 radial velocity
PARAM.volume = 0;   PARAM.volTol = 1e-14;           % impose volume conservation with first mode

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

%% SPACE DISCRETIZATION
PARAM.legendre = 1;                                     % choose 0 Chebyshev, 1 Legendre, 2 is Legendre-Lobatto
PARAM.dealiasing = 50;  ratioGrid = 2;                  % dealiasing
PARAM.tresholdCoeffs = 1e-2;     PARAM.Rescale = 1;     % option for accuracy
PARAM.n = round(ratioGrid*PARAM.dealiasing)+1;          % number of grid points

%% NEWTON METHOD PARAMETERS
PARAM.plotRes = 1;                      % plot residuals
PARAM.ResBreak = 0;                     % break when residuals are smaller
PARAM.stop = 8;                         % stop after maximum iteration
PARAM.CutStep = 1;                      % advance only partially
PARAM.plotCurv = 1;                     % plot curvature
PARAM.dh = 1e-5;                        % finite diffrence for Jacobian computation

%% CHOOSE FILENAME
PARAM.algorithm = 3;    % newton method
PARAM.filename = chooseFilename(PARAM,0,0,0);
PARAM.filenameIte = PARAM.filename;
PARAM.res = results;
PARAM.this_folder = here;
PARAM.V0 = V0;

%compute integration weight and differentiation matrices
if PARAM.legendre==1
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreIntDiff(0,1,PARAM.n+1);
    [PARAM.tcheb,PARAM.D1cheb,PARAM.D2cheb,PARAM.WGcheb] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==0
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.TTT] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==2
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreLobattoIntDiff(0,1,PARAM.n+1);
end
  
%% RUN NEWTON METHOD 
tutorial_newtonMethodSpectralXYmodesFunction(PARAM);

%the end
disp('The End')















