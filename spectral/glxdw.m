function [x,w,D] = glxdw(N)

% glxdw.m - Computation of Gauss Legendre points, weights
%           and Spectral differentiation matrix.

%  Computation of Gauss Legendre points & weights
    for N = 1:N;
     beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
     T = diag(beta,1) + diag(beta,-1); [V,D] = eig(T);
     xxx = diag(D);                       %  <- Gauss nodes
     www = 2*V(1,:).^2;                   %  <- Gauss weights
    [~,index2]=sort(xxx);
     x=xxx(index2(:));
     w=www(index2(:));
    end

    index = (1:N)';
% Construct differentiation matrix (see Fornberg book, p. 51):
  D = zeros(N,N); a = zeros(N,1);
  for k = 1:N
    notk = find(index~=k);
    a(k) = prod(x(k)-x(notk));
  end
  for k = 1:N
    notk = find(index~=k);
    D(notk,k) = (a(notk)/a(k))./(x(notk)-x(k));
    D(k,k) = sum(1./(x(k)-x(notk)));
  end
