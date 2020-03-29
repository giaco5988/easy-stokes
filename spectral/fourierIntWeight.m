%trapezi weight for Fourier integration on a priodi interval

function w = fourierIntWeight(n)

space = linspace(0,2*pi,n+1);
dl = diff(space);
dl = [dl(end) dl];

w = (dl(1:end-1)+dl(2:end))/2;