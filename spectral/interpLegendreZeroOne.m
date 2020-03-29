%legendre interpolation

function x = interpLegendreZeroOne(t,xMode)

n = numel(xMode);

PPP = lepolym(n-1,2*t-1);

manyMode = repmat(xMode,1,numel(t));

x = sum(PPP.*manyMode)';