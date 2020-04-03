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
PARAM.Ca = 0.1;                                     % capillary number
PARAM.visc = 5;     PARAM.Qdeflation = 0;           % viscosity ratio and deflation
PARAM.V0 = 4/3*pi;                                  % drop volume
PARAM.placeShapeXCM = 0;                            % when doing edge tracking, replace always droplet in x0=0 at t0
PARAM.CaNL = 0;                                     % non linear extensional flow

%% BOUNDARY CONDITIONS
PARAM.BC = 1;                                       % BC: 1 is linear extensional flow, 2 is rising droplet, 3 is nonlinear extensional flow

%% UPLOAD SHAPE FROM NEWTON METHOD RESULT
PARAM.uploadShape = 2;                  % 0 the initial shape is analytical, 1 is from continuation, 2 the initial shape is from newton Method,  3 is from edge tracking, 4 is from DNS, 5 is from Minimal Seed
PARAM.D = 0;                            % deformation parameter of the initial condition
PARAM.CaUpload = PARAM.Ca;              % capillary number of uploaded result
PARAM.viscUpload = PARAM.visc;          % viscosity ratio of uploaded result
PARAM.BCupload = PARAM.BC;              % BC of the uploaded shape
PARAM.legendreUpload = 1;               % discretization of uploaded result
PARAM.elemUpload = 50;                  % number of modes of uploaded result

%% FRAME OF REFERENCE
PARAM.dropFrame = 1;                                % choose frame of reference: 0 is lab frame, 1 is drop frame
PARAM.Unormal = 1;                                  % 0 advect nodes with lagrangian velocity 1 with normal velocity, 3 use maesh stab, 4 radial velocity
PARAM.volume = 0;   PARAM.volTol = 1e-14;           % impose volume conservation with first mode

%% SPACE DISCRETIZATION
PARAM.legendre = PARAM.legendreUpload;                  % choose 0 Chebyshev, 1 Legendre, 2 is Legendre-Lobatto
PARAM.dealiasing = PARAM.elemUpload;  ratioGrid = 2;    % dealiasing
PARAM.tresholdCoeffs = 1e-2;     PARAM.Rescale = 1;     % option for accuracy
PARAM.n = round(ratioGrid*PARAM.dealiasing)+1;          % number of grid points
PARAM.chooseMapping = 1;  PARAM.howStrong = 0.85;       % 1 is normal mapping (like spectral parametrization), 2 is clustered toward the middle (how strong set how stronly to cluster)
PARAM.chooseRemeshMapping = chooseSpectralMapping(PARAM.chooseMapping, PARAM.howStrong);

%% SET FOLDER NAMES
PARAM.res = results;
PARAM.this_folder = here;

%% compute integration weight and differentiation matrices
if PARAM.legendre==1
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreIntDiff(0,1,PARAM.n+1);
    [PARAM.tcheb,PARAM.D1cheb,PARAM.D2cheb,PARAM.WGcheb] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==0
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.TTT] = ChebyshevIntDiff(0,1,PARAM.n+1);
elseif PARAM.legendre==2
    [PARAM.t,PARAM.D1,PARAM.D2,PARAM.WG,PARAM.manyWG,PARAM.PPP] = LegendreLobattoIntDiff(0,1,PARAM.n+1);
end
  
%% PERFORM STABILITY ANALYSIS
tutorial_stabilitySpectralXYmodesFunction(PARAM)

%the end
disp('The End')















