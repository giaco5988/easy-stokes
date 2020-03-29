%compute fast fourier transorm on 2pi interval

function [x,y] = myInvFFT(theta,Xn,Yn)

    modes = numel(Yn);
    x = Xn(1);
    y = 0;
        
    for i = 1:modes

        x = x + Xn(i+1)*cos(i*theta);
        y = y + Yn(i)*sin(i*theta);

    end

end