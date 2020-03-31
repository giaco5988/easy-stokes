%given a certain value of prolate, find the corresponding delta of the
%modes (with volume always equal to a unit sphere)

function [delta,deltaF,delta_f] = find_delta(D,F)

    x = nthroot((1-D)/(1+D),3)*(1+D)/(1-D);
    y = nthroot((1-D)/(1+D),3);
    theta = 0:0.001:pi;
    
    %number of modes
    nm = 1000;
    delta_f = zeros(nm,1);
    
    shape = 1;
    
    for i = 1:nm
    
        %choose Legendre Polynomia for this mode
        P = legendre(i+1,cos(theta));
        P = P(1,:);
        
        %compute integral with trapezi rule
        int = ((sqrt((x*cos(theta(1:end-1))).^2+(y*sin(theta(1:end-1))).^2)-1).*P(1:end-1).*sin(theta(1:end-1))+...
            (sqrt((x*cos(theta(2:end))).^2+(y*sin(theta(2:end))).^2)-1).*P(2:end).*sin(theta(2:end))).*diff(theta)/2;
        
        int = sum(int);
        
%         int = @(theta) (sqrt((x*cos(theta)).^2+(y*sin(theta)).^2)-1)*P.*sin(theta);
%        delta_f(i) = 0.5*integral(int,0,pi);

        delta_f(i) = 0.5*int;
        
        shape = shape + (2*(i+1)+1)*P*delta_f(i);
    
    end
    
%     figure
%     plot(theta,shape)

    delta = norm(delta_f);
    %CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deltaF = norm(F(end-nm+1:end,1:nm)*delta_f);
    
end