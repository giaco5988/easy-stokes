%compute fast fourier transorm on 2pi interval

function [Xn,Yn] = myFFT(theta,x,modes)

    m = numel(x);
    Xn = zeros(modes+1,1);
    Yn = zeros(modes,1);

    %average
    Xn(1) = sum(x)/m;
        
    for i = 1:modes

        Xn(i+1) = 2/m*sum(x.*cos(i*theta));
        Yn(i) = 2/m*sum(x.*sin(i*theta));

    end

end