%compute fast fourier transorm on 2pi interval

function Yn = myFFTy(theta,y,modes)

    %duplicate in a smart way
    theta = [theta theta(2:end-1)+pi];
    y = [y -flip(y(2:end-1))];

    m = numel(y);
    Yn = zeros(modes,1);
        
    for i = 1:modes

        Yn(i) = 2/m*sum(y.*sin(i*theta));

    end

end