%compute Arc Length Spectrally

function l = computeTotalArcLengthSpectralCheb(x,y)

%deribvatives
xp = diff(x);  yp = diff(y);

%function to integrate
int = sqrt(xp.^2+yp.^2);

%perform integration
l = sum(int);