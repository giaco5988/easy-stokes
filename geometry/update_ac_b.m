function [ a,c ] = update_ac_b( b,lenelem,crd,ne )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
coef = lenelem.^(-1);
a = (b(2:(ne+1))-b(1:ne))./3.*coef;
c = coef.*(crd(2:(ne+1))-crd(1:ne))-1./3*lenelem.*(b(2:(ne+1))+2*b(1:ne));
end

