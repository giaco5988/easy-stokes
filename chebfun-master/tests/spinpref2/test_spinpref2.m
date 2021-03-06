% Test file for SPINPREF2/SPINPREF2:

function pass = test_spinpref2()

% Construction from STRING for SH2 equation:
pref = spinpref2('sh2');
N = pref.N;
pass(1) = isequal(N, 64);

% Construction from inputs:
pref = spinpref2('Clim', [0 10], 'dataToPlot', 'abs', 'dealias', 'off');
pass(2) = isequal(pref.Clim, [0 10]);
pass(3) = strcmpi(pref.dataToPlot, 'abs');
pass(4) = strcmpi(pref.dealias, 'off');

pref = spinpref2('dt', 1e-1, 'dtmin', 1e-10, 'dtmax', 1, 'errTol', 1e-10);
pass(5) = isequal(pref.dt, 1e-1);
pass(6) = isequal(pref.dtmin, 1e-10);
pass(7) = isequal(pref.dtmax, 1);
pass(8) = isequal(pref.errTol, 1e-10);

pref = spinpref2('iterPlot', 10, 'M', 100, 'N', 1024, 'Nmin', 10);
pass(9) = isequal(pref.iterPlot, 10);
pass(10) = isequal(pref.M, 100);
pass(11) = isequal(pref.N, 1024);
pass(12) = isequal(pref.Nmin, 10);

pref = spinpref2('Nmax', 10, 'Nplot', 2, 'plot', 'movie', 'scheme', 'lawson4');
pass(13) = isequal(pref.Nmax, 10);
pass(14) = isequal(pref.Nplot, 2);
pass(15) = strcmpi(pref.plot, 'movie');
pass(16) = strcmpi(pref.scheme, 'lawson4');

end
