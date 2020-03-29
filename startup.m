%modify default options

CURRENT_LOCATION = '~/Documents/MATLAB/';

set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

%addpath for often used directories
addpath([CURRENT_LOCATION 'physics/'])
addpath([CURRENT_LOCATION 'green_function/'])
addpath([CURRENT_LOCATION 'geometry/'])
addpath([CURRENT_LOCATION 'BEM/'])
addpath([CURRENT_LOCATION 'remesh/'])
addpath([CURRENT_LOCATION 'NonLinearFun/'])
addpath([CURRENT_LOCATION 'mex_files/'])
addpath([CURRENT_LOCATION 'spectral/'])
addpath([CURRENT_LOCATION 'finiteDifferences/'])
addpath([CURRENT_LOCATION 'chebfun-master/'])
addpath([CURRENT_LOCATION 'timeStepping/'])
addpath([CURRENT_LOCATION 'outputAndEvent/'])
addpath([CURRENT_LOCATION 'edgeTracking/'])
%addpath('splineDifferentiation/')

%Slepian library
%cd ~/Documents/MATLAB/Slepian
%initialize
%cd ~/Documents/MATLAB

%addpath for post-processing
addpath([CURRENT_LOCATION 'Postprocessing/'])
addpath([CURRENT_LOCATION 'Postprocessing/arrow'])
addpath([CURRENT_LOCATION 'Postprocessing/rising_droplet/'])
addpath([CURRENT_LOCATION 'Postprocessing/relaxation/'])
addpath([CURRENT_LOCATION 'Postprocessing/motor'])
addpath([CURRENT_LOCATION 'Postprocessing/channel_lac'])
addpath([CURRENT_LOCATION 'Postprocessing/channel_2drops/'])
addpath([CURRENT_LOCATION 'Postprocessing/channel_2D'])
addpath([CURRENT_LOCATION 'Postprocessing/extensional_flow'])
%addpath([CURRENT_LOCATION  'Postprocessing/ojwoodford-export_fig-165dc92'])
addpath([CURRENT_LOCATION 'Postprocessing/oneDropletSpectral/'])
addpath([CURRENT_LOCATION 'Postprocessing/oneDropletBEM/'])
addpath([CURRENT_LOCATION 'Postprocessing/phoretic/'])
addpath([CURRENT_LOCATION 'Postprocessing/conicalMotorSphericalBubble/'])
addpath([CURRENT_LOCATION 'Postprocessing/conicalMotorDeformableBubble/'])
addpath([CURRENT_LOCATION 'Postprocessing/verticalTube/'])
addpath([CURRENT_LOCATION 'plotFun/'])
addpath([CURRENT_LOCATION 'Postprocessing/computeFieldBEM/'])
addpath([CURRENT_LOCATION 'Postprocessing/Colormaps/'])
addpath([CURRENT_LOCATION 'Postprocessing/ojwoodford-export_fig-165dc92/'])

clear CURRENT_LOCATION






