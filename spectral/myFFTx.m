%compute fast fourier transorm on 2pi interval

function Xn = myFFTx(theta,x,modes)

    %duplicate in a smart way
    theta = [theta theta(2:end-1)+pi];
    x = [x flip(x(2:end-1))];

    m = numel(x);
    Xn = zeros(modes+1,1);

    %average
    Xn(1) = sum(x)/m;
        
    for i = 1:modes

        Xn(i+1) = 2/m*sum(x.*cos(i*theta));

    end

end