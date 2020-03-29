%compute fast fourier transorm on 2pi interval

function x = myInvFFTx(theta,Xn)

    modes = numel(Xn)-1;
    x = Xn(1);
    y = 0;
        
    for i = 1:modes

        x = x + Xn(i+1)*cos(i*theta);
        y = y + Yn(i)*sin(i*theta);

    end

end