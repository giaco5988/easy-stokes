%build evoultion operator like (Kojima 1984)

function A = my_kojima(visc,Ca)

    %kojima matrix
    
    %component of the linear operator
    n = 2:1001;
    G = 2*(visc+1)*(2*n+1)./n./(n+1).*((4*n.^4+8*n.^3+2*n.^2-8*n-6)*visc^2+...
        (8*n.^4+16*n.^3+4*n.^2-4*n+3)*visc + (4*n.^4+8*n.^3+2*n.^2+4*n));
    H = (4*n.^4+8*n.^3+14*n.^2+10*n+9)*visc^2 + (8*n.^4+16*n.^3+34*n.^2-40*n-27)*visc ...
        + (4*n.^4+8*n.^3+20*n.^2+22*n+18);
    Q = 2/Ca*(visc+1)^2*(2*n+1).^2.*(n-1).*(n+2);
    W = 2*(visc+1)*(n+2).*(n-2).*((2*n.^2+4*n+3)*visc + (2*n.^2+4*n));

    %compose the linear operator
    A = diag(-Q./G) + diag(-H(1:end-1)./G(1:end-1),1) + diag(W(2:end)./G(2:end),-1);

end