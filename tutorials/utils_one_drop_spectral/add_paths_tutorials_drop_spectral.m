function add_paths_tutorials_drop_spectral(CURRENT_LOCATION)

    set(0,'defaultaxesfontsize',25,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');

    %addpaths
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
    addpath([CURRENT_LOCATION 'Postprocessing/phoretic/'])
    addpath([CURRENT_LOCATION 'Postprocessing'])
    addpath([CURRENT_LOCATION 'Postprocessing/computeFieldBEM'])
    addpath([CURRENT_LOCATION 'plotFun'])
    
end