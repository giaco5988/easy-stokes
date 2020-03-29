%legendre interpolation

function x = interpLegendre(t,x,In,End,Mode)

L = End-In;
n = numel(Mode);

PPP = lepolym(n,t);