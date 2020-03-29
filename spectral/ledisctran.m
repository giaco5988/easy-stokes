% The function s=ledisctran(n,x,p,iflag) performs the discrete Legendre transforms 
% between the physical space (i.e., physical values) and frequence space
% (Legendre expansion coefficients) at the Legendre-Gauss-Lobatto points
% Input:
%  n, x,w--- number of LGL points in x, where (x,w) can be computed by
%          [x,w]=legslb(n). Note: x,w are column vectors  
%  iflag==0--- forward transform  
%    p--- (input) physical values at collocation points
%    s--- (output) expansion coefficients 
%  iflag not= 0--- backward transform  
%    p--- (input) expansion coefficients 
%    s--- (output) physical values at collocation points 
%
%  See Page 101 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepolym()   
%  Last modified on August 31, 2011


function s=ledisctran(n,x,w,p,iflag)
 T=lepolym(n-1,x); % compute the Legndre polynomials up to order n-1. Note: T(i,j)=L_i(x_j) 
 if iflag==0, 
     s=(T*(p.*w)).*[[0:n-2]'+0.5;(n-1)/2]; % see (3.193)
     return;
 end
 
 s=T'*p;  % see (3.194)
 return
 
 

